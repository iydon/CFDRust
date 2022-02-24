! Using the department Unix computer and a simple Fortran code, investigate how an input real number is outputted (i.e., interpreted by a computer in bits and bytes) under the single precision representation.
!
! (a) Input = 0.61E-45
! (b) Input = 0.82E-45
! (c) Input = 2.13E-45
! (d) Input = 2.34E-45
! (e) Input = 3.65E-10
! (f) Input = 4.06E+30
! (g) Input = 3.40E+38
! (h) Input = 3.41E+38
!
! Namely, in each case, (i) report the actual output of your input; (ii) explain the precise binary form used by the computer. Your code is expected to be less than 10 lines.
!
! Hint: you may consult the Wiki page at https://en.wikipedia.org/wiki/Single-precision_floating-point_format

program main
   use, intrinsic :: iso_fortran_env, only: sp=>real32
   implicit none

   real(sp) :: number

   write(*, fmt='(A)', advance='no') 'Input  = '
   read(*, *) number
   write(*, fmt='(A, ES14.8, A, b32.32)') 'Output = ', number, ' = ', transfer(number, 0)
end program main

! #[compile]
! {
!     "fortran": "gfortran",
!     "options": []
! }

! #[test]
! [
!     {
!         "stdin": "0.61E-45",
!         "stdout": "Input  = Output = 0.00000000E+00 = 00000000000000000000000000000000\n"
!     },
!     {
!         "stdin": "0.82E-45",
!         "stdout": "Input  = Output = 1.40129846E-45 = 00000000000000000000000000000001\n"
!     },
!     {
!         "stdin": "2.13E-45",
!         "stdout": "Input  = Output = 2.80259693E-45 = 00000000000000000000000000000010\n"
!     },
!     {
!         "stdin": "2.34E-45",
!         "stdout": "Input  = Output = 2.80259693E-45 = 00000000000000000000000000000010\n"
!     },
!     {
!         "stdin": "3.65E-10",
!         "stdout": "Input  = Output = 3.65000002E-10 = 00101111110010001010100100101111\n"
!     },
!     {
!         "stdin": "4.06E+30",
!         "stdout": "Input  = Output = 4.05999996E+30 = 01110010010011001111101001000101\n"
!     },
!     {
!         "stdin": "3.40E+38",
!         "stdout": "Input  = Output = 3.39999995E+38 = 01111111011111111100100110011110\n"
!     },
!     {
!         "stdin": "3.41E+38",
!         "stdout": "Input  = Output =       Infinity = 01111111100000000000000000000000\n"
!     }
! ]
