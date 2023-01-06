program cholesky_mpi_block
   use mpi_f08
   implicit none
   integer, parameter :: m = 10
   integer :: i, j, k, n, nb, t, s, p
   real*8, allocatable :: id(:,:), a(:,:), l(:,:), l_21(:,:), l_massage(:)
   real*8, allocatable :: b(:,:), x(:,:), y(:,:)
   real*8 :: time  ! time counter
   integer :: nproc, rank, ierr
   type(MPI_Status) :: status
   ! Interface for subroutines
   interface
      subroutine cholesky_serial(n, l)
         implicit none
         integer, intent(in) :: n
         real*8 :: l(n,n)
         real*8 :: s
         integer :: i, ip, j
      end subroutine cholesky_serial

      subroutine updating_l_21(nrows, ncols, l_21_, l_11_)
         integer, intent(in) :: nrows, ncols
         integer :: i_, j_, k_
         real*8 :: l_21_(nrows, ncols), x(nrows, ncols)
         real*8, intent(in) :: l_11_(ncols, ncols)
      end subroutine
   end interface
   ! Initialization of the MPI environment
   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
   if (nproc .ne. 1) then
      write(*,*) 'nproc should be 1'
      call exit(1)
   endif
   ! Reading n and calculation nb
   if (rank .eq. 0) then
      read(*, *) n, nb
   end if
   ! Initialization of the time counter
   if (rank .eq. 0) then
      time = MPI_WTIME()
   end if
   ! Id matrix
   allocate(id(n, n))
   ForAll(i = 1:n, j = 1:n) id(i,j) = (i/j)*(j/i)
   ! Initialization of the matrix A
   allocate(a(n, n))
   a(:, :) = 1.0
   a = (a + transpose(a)) / 2 + n * id
   ! Initialization of the matrix L and distributing it
   allocate(l(n,n))
   l(:,:) = 0.0
   do i = 1, n, nb
      j = min(i + nb - 1, n)
      l(i:n, i:j) = a(i:n, i:j)
   end do
   ! Block version of the Cholesky decomposition
   do i = 1, n, nb
      j = min(i + nb - 1, n)
      ! Updating the (i:n, i:j) block of the L matrix
      if (mod( (i-1)/nb, nproc ) .eq. rank) then
         ! Updating the (i:j, i:j) block of the L matrix
         call cholesky_serial(j - i + 1, l(i:j, i:j))
      end if
      if (j .lt. n) then
         if (mod( (i-1)/nb, nproc ) .eq. rank) then
            ! Updating the (j+1:n, i:j) block of the matrix L
            call updating_l_21(n - j, j - i + 1, l(j + 1:n, i:j), l(i:j, i:j))
         end if
         ! Allocation the useful array
         allocate(l_21(j+1:n, i:j))
         l_21 = l(j + 1:n, i:j)
         ! Updating the other blocks of the L matrix
         do k = i + nb, n, nb
            p = min(k + nb - 1, n)
            ! Updating the (k:n, k:p) block
            l(k:n, k:p) = l(k:n, k:p) - matmul(l_21(k:n, i:j), transpose(l_21(k:p, i:j)))
         end do
         deallocate(l_21)
      end if
   end do
   ! Showing the transposed matrix L
   ! do i = rank * nb + 1, n, nb * nproc
   !    j = min(i + nb - 1, n)
   !    do k = i, j
   !       print *, k, l(:, k)
   !    end do
   ! end do
   ! Initialization of matrices B, X and Y
   allocate(b(n,m), x(n,m), y(n,m))
   do i = 1, m
      b(:,i) = i
   end do
   ! Backward substitution
   do i = 1, m
      ! Finding y(:,i)
      y(1,i) = b(1,i) / l(1,1)
      do j = 2, n
         y(j,i) = b(j,i)
         do k = 1, j - 1
            y(j,i) = y(j,i) - l(j,k) * y(k,i)
         end do
         y(j,i) = y(j,i) / l(j,j)
      end do
      ! Finding x(:,i)
      x(n,i) = y(n,i) / l(n,n)
      do j = n-1, 1, -1
         x(j,i) = y(j,i)
         do k = j+1, n
            x(j,i) = x(j,i) - l(k,j) * x(k,i)
         end do
         x(j,i) = x(j,i) / l(j,j)
      end do 
   end do
   ! Showing the transposed matrix X
   ! do i = 1, m
   !    if (mod(i-1, nproc) .eq. rank) then
   !       print *, i, x(:,i)
   !    end if
   ! end do
   ! Time measurement
   if (rank .eq. 0) then
      time = MPI_WTIME() - time
      print *, n, time
   end if
   deallocate(id, a, l)
   deallocate(b, x, y)
   call MPI_FINALIZE(ierr)
end program cholesky_mpi_block

subroutine cholesky_serial(n, l)
   implicit none
   integer, intent(in) :: n
   real*8 :: l(n,n)
   real*8 :: s
   integer :: i, ip, j
   ! Serial Cholesky decomposition
   do  i = 1, n
      s = l(i,i)
      do  ip = 1, i - 1
         s = s - l(i,ip) * l(i,ip)
      end do
      l(i,i) = sqrt(s)
      do  j = i + 1, n
         s = l(j,i)
         do  ip = 1, i-1
            s = s - l(i,ip) * l(j,ip)
         end do
         l(j,i) = s / l(i,i)
      end do
   end do
   do i = 1, n
      l(i,i+1:n) = 0.0
   end do
end subroutine cholesky_serial

subroutine updating_l_21(nrows, ncols, l_21_, l_11_)
   integer, intent(in) :: nrows, ncols
   integer :: i_, j_, k_
   real*8 :: l_21_(nrows, ncols), x(nrows, ncols)
   real*8, intent(in) :: l_11_(ncols, ncols)
   do i_ = 1, nrows
      x(i_, 1) = l_21_(i_, 1) / l_11_(1, 1)
      do j_ = 2, ncols
         x(i_, j_) = l_21_(i_, j_)
         do k_ = 1, j_ - 1
            x(i_, j_) = x(i_, j_) - x(i_, k_) * l_11_(j_, k_)
         end do
         x(i_, j_) = x(i_, j_) / l_11_(j_, j_)
      end do
   end do
   l_21_ = x
end subroutine