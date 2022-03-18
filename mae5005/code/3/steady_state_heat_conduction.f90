program main
   use, intrinsic :: iso_fortran_env, only: f64=>real64
   implicit none
   ! user-specified variables
   integer :: nx, ny
   ! variable
   integer :: ith, jth
   real(f64) :: x, y
   real(f64), allocatable :: field(:, :)
   ! function
   real(f64) :: T

   write(*, fmt='(A)', advance='no') 'nx = '; read(*, *) nx
   write(*, fmt='(A)', advance='no') 'ny = '; read(*, *) ny
   allocate(field(nx+1, ny+1))
   open(1, file='numerical.dat', action='write')

   do ith = 0, nx
      x = 1.0 * ith / nx
      do jth = 0, ny
         y = 1.0 * jth / ny
         field(ith+1, jth+1) = T(x, y, 999)
      end do
   end do
   write(1, *) field

   close(1)
   deallocate(field)
end program main

function T(x, y, number) result(ans)
   use, intrinsic :: iso_fortran_env, only: f64=>real64
   ! input and output
   real(f64), intent(in) :: x, y
   integer, intent(in) :: number
   real(f64) :: ans
   ! function
   real(f64) :: term
   ! variable
   integer :: n
   real(f64) :: delta

   ans = 0.0
   do n = 1, number
      delta = term(x, y, n)
      if (isnan(delta)) then
         write(*, fmt='(I3)') n-1
         exit
      end if
      ans = ans + delta
   end do
end function T

function term(x, y, n) result(ans)
   use, intrinsic :: iso_fortran_env, only: f64=>real64
   ! input and output
   real(f64), intent(in) :: x, y
   integer, intent(in) :: n
   real(f64) :: ans
   ! parameter
   real(f64), parameter :: PI = 4.0_f64*datan(1.0_f64)
   real(f64) :: m

   m = (2.0*n - 1.0) * PI
   ans = (4.0 * dsin(m*x)) &
      * (dexp(-2*m*y) - dexp(-2*m)) &
      / (m * dexp(-m*y) * (1.0-dexp(-2*m)))
end function term

! #[compile]
! {
!     "fortran": "gfortran",
!     "options": []
! }

! #[test]
! [
!     {
!         "stdin": "9\n9",
!         "stdout": "nx = ny = 534\n356\n267\n213\n178\n152\n133\n119\n534\n356\n267\n213\n178\n152\n133\n119\n534\n356\n267\n213\n178\n152\n133\n119\n534\n356\n267\n213\n178\n152\n133\n119\n534\n356\n267\n213\n178\n152\n133\n119\n534\n356\n267\n213\n178\n152\n133\n119\n534\n356\n267\n213\n178\n152\n133\n119\n534\n356\n267\n213\n178\n152\n133\n119\n534\n356\n267\n213\n178\n152\n133\n119\n534\n356\n267\n213\n178\n152\n133\n119\n"
!     }
! ]
