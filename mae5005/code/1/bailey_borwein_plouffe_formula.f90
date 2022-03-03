program main
   use, intrinsic :: iso_fortran_env, only: f32=>real32, f64=>real64, f128=>real128
   implicit none

   CHARACTER(LEN=23) :: format = '(A13, I4, A7, f36.34)'
   integer :: n32, n64, n128
   real(f32) :: pi32
   real(f64) :: pi64
   real(f128) :: pi128

   call calc_pi_f32(n32, pi32)
   write(*, format) '32-bit: n = ', n32, ', pi = ', pi32

   call calc_pi_f64(n64, pi64)
   write(*, format) '64-bit: n = ', n64, ', pi = ', pi64

   call calc_pi_f128(n128, pi128)
   write(*, format) '128-bit: n = ', n128, ', pi = ', pi128

   write(*, '(A)') '                   pi = 3.1415926535897932384626433832795029'
end program main

subroutine calc_pi_f32(n, pi)
   use, intrinsic :: iso_fortran_env, only: f32=>real32

   integer, intent(out) :: n
   real(f32), intent(out) :: pi

   real(f32) :: old, a, b

   n = 1 ; pi = 0.0 ; old = pi ; a = 1.0 ; b = 0.0
   do while (.true.)
      pi = pi + (4.0/(b+1.0) - 2.0/(b+4.0) - 1.0/(b+5.0) - 1.0/(b+6.0)) / a
      if (old == pi) then
         exit
      end if
      n = n + 1 ; old = pi ; a = a * 16.0 ; b = b + 8.0
   end do
end subroutine calc_pi_f32

subroutine calc_pi_f64(n, pi)
   use, intrinsic :: iso_fortran_env, only: f64=>real64

   integer, intent(out) :: n
   real(f64), intent(out) :: pi

   real(f64) :: old, a, b

   n = 1 ; pi = 0.0 ; old = pi ; a = 1.0 ; b = 0.0
   do while (.true.)
      pi = pi + (4.0/(b+1.0) - 2.0/(b+4.0) - 1.0/(b+5.0) - 1.0/(b+6.0)) / a
      if (old == pi) then
         exit
      end if
      n = n + 1 ; old = pi ; a = a * 16.0 ; b = b + 8.0
   end do
end subroutine calc_pi_f64

subroutine calc_pi_f128(n, pi)
   use, intrinsic :: iso_fortran_env, only: f128=>real128

   integer, intent(out) :: n
   real(f128), intent(out) :: pi

   real(f128) :: old, a, b

   n = 1 ; pi = 0.0 ; old = pi ; a = 1.0 ; b = 0.0
   do while (.true.)
      pi = pi + (4.0/(b+1.0) - 2.0/(b+4.0) - 1.0/(b+5.0) - 1.0/(b+6.0)) / a
      if (old == pi) then
         exit
      end if
      n = n + 1 ; old = pi ; a = a * 16.0 ; b = b + 8.0
   end do
end subroutine calc_pi_f128

! #[compile]
! {
!     "fortran": "gfortran",
!     "options": []
! }

! #[test]
! [
!     {
!         "stdout": " 32-bit: n =    6, pi = 3.1415925025939941406250000000000000\n 64-bit: n =   12, pi = 3.1415926535897931159979634685441852\n128-bit: n =   27, pi = 3.1415926535897932384626433832795024\n                   pi = 3.1415926535897932384626433832795029\n"
!     }
! ]
