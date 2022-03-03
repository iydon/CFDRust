module solution
   use, intrinsic :: iso_fortran_env, only: f64=>real64
   implicit none

   public stepwise, numerical, analytical
contains
   subroutine stepwise(N, L, nu, u0, dt, u)
      integer, intent(in) :: N
      real(f64), intent(in) :: L, nu, u0, dt
      real(f64), dimension(N+1), intent(inout) :: u

      u(2: N) = u(2: N) + dt * ( &
         + 2.0*nu*u0/L**2 & ! u₀ = -L²/(2ρν) pₓ
         + nu/(2.0*L/N)**2 * (u(1: N-1)-2.0*u(2: N)+u(3: N+1)) &
      )
   end subroutine stepwise

   subroutine numerical(N, nt, L, nu, u0, dt, u)
      integer, intent(in) :: N, nt
      real(f64), intent(in) :: L, nu, u0, dt
      real(f64), dimension(N+1), intent(out) :: u

      integer :: ith

      u = 0.0 ! initial condition
      do ith = 1, nt
         call stepwise(N, L, nu, u0, dt, u)
      end do
   end subroutine numerical

   subroutine analytical(N, nt, L, nu, u0, dt, u, term)
      integer, intent(in) :: N, nt
      real(f64), intent(in) :: L, nu, u0, dt
      integer, intent(in) :: term ! when term<0, add until it stops increasing
      real(f64), dimension(N+1), intent(out) :: u

      real(f64), parameter :: PI = 4.0 * atan(1.0)
      integer :: ith, k
      real(f64) :: t, y, dy, C, old

      t = nt * dt; y = -L; dy = 2.0 * L / N
      do ith = 1, N+1
         u(ith) = 1.0 - y**2 / L**2
         if (term >= 0) then
            do k = 0, term-1
               C = PI * (k + 0.5)
               u(ith) = u(ith) - 4.0*(-1.0)**k/C**3 * exp(-C**2*nu*t/L**2) * cos(C*y/L)
            end do
         else
            k = 0
            do while (.true.)
               C = PI * (k + 0.5); old = u(ith)
               u(ith) = u(ith) - 4.0*(-1.0)**k/C**3 * exp(-C**2*nu*t/L**2) * cos(C*y/L)
               if (u(ith) == old) then; exit; end if
               k = k + 1
            end do
         end if
         u(ith) = u(ith) * u0; y = y + dy
      end do
   end subroutine analytical
end module solution

program main
   use, intrinsic :: iso_fortran_env, only: f64=>real64
   use solution, only: stepwise, numerical, analytical
   implicit none

   ! user-specified variables
   integer :: N ! N = 8, 16, 32
   real(f64) :: nutL ! ν t / L² = 0.2, 1.0, 10.0
   ! known parameters
   real(f64), parameter :: L = 1.0, nu = 0.1
   real(f64), parameter :: u0 = 1.0 ! u₀ = -L² / (2 ρ ν) pₓ
   ! derived variables
   integer :: nt ! ν nt Δt / L² = nutL
   real(f64) :: dy ! Δy = 2 L / N
   real(f64) :: dt ! Δt = 0.32 Δy² / ν
   ! numerical and analytical solutions
   real(f64), allocatable :: un(:), ua(:)
   ! unimportant temporary variables
   integer :: ith, term = 100
   character(len=25) :: format = '(A, I4, A, f9.6, A, f8.6)'

   write(*, fmt='(A)', advance='no') 'N     = '; read(*, *) N
   write(*, fmt='(A)', advance='no') 'νt/L² = '; read(*, *) nutL
   dy = 2.0 * L / N
   dt = 0.32 * dy**2 / nu ! 1.28 * L**2 / N**2 / nu
   nt = int(nutL * L**2 / nu / dt) ! int(nutL / 1.28 * N**2)
   allocate(un(N+1), ua(N+1))
   open(1, file='numerical.dat', action='write')
   open(2, file='analytical.dat', action='write')

   ! plot and compare the velocity profiles from different grid resolutions
   call numerical(N, nt, L, nu, u0, dt, un)
   call analytical(N, nt, L, nu, u0, dt, ua, term)
   write(1, *) un
   write(2, *) ua
   ! find out at what dimensionless time, νt/L², the velocity at the center of the channel (y=0) reaches the value of 0.99u₀
   write(*, format) '[C] nt, νt/L², u(0, T)/u₀ = ', nt, ', ', nutL, ', ', un(int(N/2)+1)/u0
   ith = 0; un = 0.0
   do while (.true.)
      ith = ith + 1
      call stepwise(N, L, nu, u0, dt, un)
      if (un(int(N/2)+1) >= 0.99*u0) then
         write(*, format) '[N] nt, νt/L², u(0, T)/u₀ = ', ith, ', ', nu*ith*dt/L**2, ', ', un(int(N/2)+1)/u0
         exit
      end if
   end do
   ith = 0
   do while (.true.)
      ith = ith + 1
      call analytical(N, ith, L, nu, u0, dt, ua, term)
      if (ua(int(N/2)+1) >= 0.99*u0) then
         write(*, format) '[A] nt, νt/L², u(0, T)/u₀ = ', ith, ', ', nu*ith*dt/L**2, ', ', ua(int(N/2)+1)/u0
         exit
      end if
   end do

   close(1); close(2)
   deallocate(un, ua)
end program

! #[compile]
! {
!     "fortran": "gfortran",
!     "options": []
! }

! #[test]
! [
!     {
!         "stdin": "8\n0.2",
!         "stdout": "N     = νt/L² = [C] nt, νt/L², u(0, T)/u₀ =   10,  0.200000, 0.373990\n[N] nt, νt/L², u(0, T)/u₀ =   93,  1.860000, 0.990081\n[A] nt, νt/L², u(0, T)/u₀ =   94,  1.880000, 0.990020\n"
!     },
!     {
!         "stdin": "8\n1.0",
!         "stdout": "N     = νt/L² = [C] nt, νt/L², u(0, T)/u₀ =   50,  1.000000, 0.915054\n[N] nt, νt/L², u(0, T)/u₀ =   93,  1.860000, 0.990081\n[A] nt, νt/L², u(0, T)/u₀ =   94,  1.880000, 0.990020\n"
!     },
!     {
!         "stdin": "8\n10.0",
!         "stdout": "N     = νt/L² = [C] nt, νt/L², u(0, T)/u₀ =  500, 10.000000, 1.000000\n[N] nt, νt/L², u(0, T)/u₀ =   93,  1.860000, 0.990081\n[A] nt, νt/L², u(0, T)/u₀ =   94,  1.880000, 0.990020\n"
!     },
!     {
!         "stdin": "16\n0.2",
!         "stdout": "N     = νt/L² = [C] nt, νt/L², u(0, T)/u₀ =   40,  0.200000, 0.371261\n[N] nt, νt/L², u(0, T)/u₀ =  375,  1.875000, 0.990034\n[A] nt, νt/L², u(0, T)/u₀ =  376,  1.880000, 0.990020\n"
!     },
!     {
!         "stdin": "16\n1.0",
!         "stdout": "N     = νt/L² = [C] nt, νt/L², u(0, T)/u₀ =  200,  1.000000, 0.913117\n[N] nt, νt/L², u(0, T)/u₀ =  375,  1.875000, 0.990034\n[A] nt, νt/L², u(0, T)/u₀ =  376,  1.880000, 0.990020\n"
!     },
!     {
!         "stdin": "16\n10.0",
!         "stdout": "N     = νt/L² = [C] nt, νt/L², u(0, T)/u₀ = 2000, 10.000000, 1.000000\n[N] nt, νt/L², u(0, T)/u₀ =  375,  1.875000, 0.990034\n[A] nt, νt/L², u(0, T)/u₀ =  376,  1.880000, 0.990020\n"
!     },
!     {
!         "stdin": "32\n0.2",
!         "stdout": "N     = νt/L² = [C] nt, νt/L², u(0, T)/u₀ =  160,  0.200000, 0.370603\n[N] nt, νt/L², u(0, T)/u₀ = 1503,  1.878750, 0.990023\n[A] nt, νt/L², u(0, T)/u₀ = 1504,  1.880000, 0.990020\n"
!     },
!     {
!         "stdin": "32\n1.0",
!         "stdout": "N     = νt/L² = [C] nt, νt/L², u(0, T)/u₀ =  800,  1.000000, 0.912637\n[N] nt, νt/L², u(0, T)/u₀ = 1503,  1.878750, 0.990023\n[A] nt, νt/L², u(0, T)/u₀ = 1504,  1.880000, 0.990020\n"
!     },
!     {
!         "stdin": "32\n10.0",
!         "stdout": "N     = νt/L² = [C] nt, νt/L², u(0, T)/u₀ = 8000, 10.000000, 1.000000\n[N] nt, νt/L², u(0, T)/u₀ = 1503,  1.878750, 0.990023\n[A] nt, νt/L², u(0, T)/u₀ = 1504,  1.880000, 0.990020\n"
!     }
! ]
