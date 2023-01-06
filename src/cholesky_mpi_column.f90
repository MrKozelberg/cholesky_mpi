program cholesky_mpi_column
   use mpi_f08
   implicit none
   integer, parameter :: m = 10
   integer :: i, j, k, n
   real*4, allocatable :: id(:,:), a(:,:), l(:,:), l_column(:)
   real*8 :: time  ! time counter
   integer :: nproc, rank, ierr
   type(MPI_Status) :: status
   ! Initialization of the MPI environment
   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
   if (rank .eq. 0) then
      time = MPI_WTIME()
   end if
   ! Reading n
   if (rank .eq. 0) then
      read(*,*) n
   end if
   ! Broadcasting n to other processes
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   ! Id matrix
   allocate(id(n,n))
   ForAll(i = 1:n, j = 1:n) id(i,j) = (i/j)*(j/i)
   ! Initialization of the matrix A
   allocate(a(n,n))
   a(:,:) = 1.0
   a = (a + transpose(a)) / 2 + n * id
   ! ! Printing the matrix A
   ! if (rank.eq.0) then
   !    print *, "A"
   !    do i = 1, n
   !       print *, a(:,i)
   !    end do
   ! end if
   ! Initialization of the matrix L
   allocate(l(n,n))
   l(:,:) = 0.0
   do k = rank + 1, n, nproc
      l(k,k:n) = a(k,k:n)
   end do
   ! Cholesky decomposition
   do k = 1, n
      ! Allocation the vector l_column
      allocate(l_column(n-k))
      ! Update the k-th column of the L matrix
      if (mod(k-1,nproc).eq.rank) then
         l(k,k) = sqrt(l(k,k))
         l(k,k+1:n) = l(k,k+1:n) / l(k,k)
         l_column(:) = l(k,k+1:n)
      end if
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      ! Broadcasting the new column to other processes
      call MPI_BCAST(l_column, n-k, MPI_REAL, mod(k-1,nproc), MPI_COMM_WORLD, ierr)
      if (k .eq. 1) then
         print*, rank, l_column
      end if
      ! Update rest columns
      do j = k + 1, n
         if (mod(j-1,nproc).eq.rank) then
            l(j,j:n) = l(j,j:n) - l_column(j-k:n) * l_column(j-k)
         end if
      end do
      ! Deallocation the vector l_column
      deallocate(l_column)
   end do
   ! Time measurement
   if (rank .eq. 0) then
      time = MPI_WTIME() - time
      print *, n, time
      print *, 'L^T'
   end if
   ! ! Printing the result
   ! do k = 1, n
   !    if (mod(k-1,nproc).eq.rank) then
   !       print *, l(k,:)
   !    end if
   ! end do
   deallocate(id)
   deallocate(a)
   deallocate(l)
   call MPI_FINALIZE(ierr)
end program cholesky_mpi_column