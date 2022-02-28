! In 1995 (very recent in terms of the history of mathematics), Simon Plouffe discovered an amazing formulae for computing π, which, in principle, allows one to compute π to any binary (decimal) digit. His formulae is
!
! \pi = \sum_{k=0}^{\infty} \left[
!     \frac{1}{16^k} \left(
!         \frac{4}{8k+1} - \frac{2}{8k+4} - \frac{1}{8k+5} - \frac{1}{8k+6}
!     \right)
! \right]
!
! Clearly, the different terms in terms of k, decreases quickly as k is increased. Here, we shall use the above formulae to compute π, at different accuracies and investigate how many terms are needed for each case.
!
! (a) Write a single precision Fortran code, to compute π using the above formulae;
! (b) Write a double precision Fortran code, to compute π using the above formulae;
! (c) Write a quadruple precision Fortran code, to compute π using the above formulae.
!
! Report, in each case, what do you obtain? what is the accuracy? How many terms, in terms of k, are needed to achieve the accuracy? Your code for each part is expected to be less than 20 lines. The value of π, up to 100th decimal places, is given here
!
! π = 3.1415926535 8979323846 2643383279 5028841971 6939937510
!       5820974944 5923078164 0628620899 8628034825 3421170679

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
