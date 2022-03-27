module solution
   use, intrinsic :: iso_fortran_env, only: f64=>real64
   implicit none

   public stepwise, numerical, analytical
contains
   subroutine stepwise(N, L, nu, u0, dt, u)
      use, intrinsic :: iso_fortran_env, only: f64=>real64
   
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
      integer, intent(in) :: term
      real(f64), dimension(N+1), intent(out) :: u

      real(f64), parameter :: PI = 4.0 * datan(1.0_f64)
      integer :: ith, k
      real(f64) :: t, y, dy, C

      t = nt*dt; y = -L; dy = 2.0*L/N
      do ith = 1, N+1
         u(ith) = 1.0 - y**2 / L**2
         do k = 0, term-1
            C = PI * (k + 0.5)
            u(ith) = u(ith) - 4.0*(-1.0)**k/C**3 * dexp(-C**2*nu*t/L**2) * dcos(C*y/L)
         end do
         u(ith) = u(ith) * u0; y = y + dy
      end do
   end subroutine analytical

   subroutine modified(N, nt, L, nu, u0, dt, u, term)
      integer, intent(in) :: N, nt
      real(f64), intent(in) :: L, nu, u0, dt
      integer, intent(in) :: term
      real(f64), dimension(N+1), intent(out) :: u

      real(f64), parameter :: PI = 4.0 * datan(1.0_f64)
      integer :: ith, k
      real(f64) :: t, y, dy, C, alpha, num

      t = nt*dt; y = -L; dy = 2.0*L/N; alpha = nu*dt/dy**2
      do ith = 1, N+1
         u(ith) = 1.0 - y**2 / L**2
         do k = 0, term-1
            C = PI * (k + 0.5)
            num = nu * ( &
               + 1.0 &
               + (alpha/2.0-1.0/12.0)*(2.0*C/N)**2 &
               + (alpha**2/3.0-alpha/12.0+1.0/360.0)*(2.0*C/N)**4 &
            )
            u(ith) = u(ith) - 4.0*(-1.0)**k/C**3 * dexp(-C**2*num*t/L**2) * dcos(C*y/L)
         end do
         u(ith) = u(ith) * u0; y = y + dy
      end do
   end subroutine modified
end module solution

program main
   use, intrinsic :: iso_fortran_env, only: f64=>real64
   use solution, only: stepwise, numerical, analytical, modified
   implicit none

   ! user-specified variables
   integer :: N, nt
   real(f64) :: L, nu, dt
   ! known parameters
   integer, parameter :: term = 100
   real(f64), parameter :: u0 = 1.0 ! u/u₀ ≡ u
   ! derived variables
   real(f64) :: dy
   ! numerical and analytical solutions
   real(f64), allocatable :: un(:), ua(:), um(:)
   ! unimportant temporary variables
   integer :: ith
   character(len=25) :: format = '(A, f9.7, A, f9.7)'

   write(*, fmt='(A)', advance='no') 'N  = '; read(*, *) N
   write(*, fmt='(A)', advance='no') 'nt = '; read(*, *) nt
   write(*, fmt='(A)', advance='no') 'L  = '; read(*, *) L
   write(*, fmt='(A)', advance='no') 'ν  = '; read(*, *) nu
   write(*, fmt='(A)', advance='no') 'Δt = '; read(*, *) dt
   dy = 2.0 * L / N
   allocate(un(N+1), ua(N+1), um(N+1))
   open(1, file='numerical.dat', action='write')
   open(2, file='analytical.dat', action='write')
   open(3, file='modified.dat', action='write')

   ith = 0; un = 0.0
   do ith = 1, nt
      call stepwise(N, L, nu, u0, dt, un)
      call analytical(N, ith, L, nu, u0, dt, ua, term)
      call modified(N, ith, L, nu, u0, dt, um, term)
      write(1, *) un
      write(2, *) ua
      write(3, *) um
      write(*, format) 'uₘ-u, uₙ-u = ', dabs(um(N/2+1)-ua(N/2+1)), ', ', dabs(un(N/2+1)-ua(N/2+1))
   end do

   close(1); close(2); close(3)
   deallocate(un, ua, um)
end program

! #[compile]
! {
!     "fortran": "gfortran",
!     "options": []
! }

! #[test]
! [
!     {
!         "stdin": "16\n400\n1.0\n0.1\n0.049999999999999996",
!         "stdout": "N  = nt = L  = ν  = Δt = uₘ-u, uₙ-u = 0.0000001, 0.0000000\nuₘ-u, uₙ-u = 0.0000001, 0.0000000\nuₘ-u, uₙ-u = 0.0000002, 0.0000000\nuₘ-u, uₙ-u = 0.0000002, 0.0000000\nuₘ-u, uₙ-u = 0.0000001, 0.0000001\nuₘ-u, uₙ-u = 0.0000000, 0.0000005\nuₘ-u, uₙ-u = 0.0000016, 0.0000023\nuₘ-u, uₙ-u = 0.0000069, 0.0000077\nuₘ-u, uₙ-u = 0.0000169, 0.0000178\nuₘ-u, uₙ-u = 0.0000321, 0.0000330\nuₘ-u, uₙ-u = 0.0000522, 0.0000532\nuₘ-u, uₙ-u = 0.0000768, 0.0000779\nuₘ-u, uₙ-u = 0.0001051, 0.0001063\nuₘ-u, uₙ-u = 0.0001362, 0.0001377\nuₘ-u, uₙ-u = 0.0001696, 0.0001712\nuₘ-u, uₙ-u = 0.0002045, 0.0002063\nuₘ-u, uₙ-u = 0.0002403, 0.0002423\nuₘ-u, uₙ-u = 0.0002766, 0.0002788\nuₘ-u, uₙ-u = 0.0003130, 0.0003154\nuₘ-u, uₙ-u = 0.0003492, 0.0003517\nuₘ-u, uₙ-u = 0.0003848, 0.0003875\nuₘ-u, uₙ-u = 0.0004198, 0.0004226\nuₘ-u, uₙ-u = 0.0004539, 0.0004568\nuₘ-u, uₙ-u = 0.0004870, 0.0004901\nuₘ-u, uₙ-u = 0.0005191, 0.0005223\nuₘ-u, uₙ-u = 0.0005501, 0.0005534\nuₘ-u, uₙ-u = 0.0005800, 0.0005833\nuₘ-u, uₙ-u = 0.0006087, 0.0006121\nuₘ-u, uₙ-u = 0.0006363, 0.0006398\nuₘ-u, uₙ-u = 0.0006627, 0.0006663\nuₘ-u, uₙ-u = 0.0006881, 0.0006916\nuₘ-u, uₙ-u = 0.0007123, 0.0007159\nuₘ-u, uₙ-u = 0.0007355, 0.0007391\nuₘ-u, uₙ-u = 0.0007577, 0.0007613\nuₘ-u, uₙ-u = 0.0007788, 0.0007824\nuₘ-u, uₙ-u = 0.0007990, 0.0008026\nuₘ-u, uₙ-u = 0.0008182, 0.0008219\nuₘ-u, uₙ-u = 0.0008366, 0.0008402\nuₘ-u, uₙ-u = 0.0008541, 0.0008577\nuₘ-u, uₙ-u = 0.0008707, 0.0008743\nuₘ-u, uₙ-u = 0.0008866, 0.0008902\nuₘ-u, uₙ-u = 0.0009017, 0.0009053\nuₘ-u, uₙ-u = 0.0009160, 0.0009196\nuₘ-u, uₙ-u = 0.0009297, 0.0009332\nuₘ-u, uₙ-u = 0.0009427, 0.0009462\nuₘ-u, uₙ-u = 0.0009550, 0.0009585\nuₘ-u, uₙ-u = 0.0009668, 0.0009702\nuₘ-u, uₙ-u = 0.0009779, 0.0009813\nuₘ-u, uₙ-u = 0.0009885, 0.0009919\nuₘ-u, uₙ-u = 0.0009985, 0.0010019\nuₘ-u, uₙ-u = 0.0010080, 0.0010113\nuₘ-u, uₙ-u = 0.0010170, 0.0010203\nuₘ-u, uₙ-u = 0.0010255, 0.0010288\nuₘ-u, uₙ-u = 0.0010336, 0.0010368\nuₘ-u, uₙ-u = 0.0010412, 0.0010444\nuₘ-u, uₙ-u = 0.0010484, 0.0010516\nuₘ-u, uₙ-u = 0.0010552, 0.0010583\nuₘ-u, uₙ-u = 0.0010615, 0.0010647\nuₘ-u, uₙ-u = 0.0010675, 0.0010706\nuₘ-u, uₙ-u = 0.0010732, 0.0010762\nuₘ-u, uₙ-u = 0.0010785, 0.0010815\nuₘ-u, uₙ-u = 0.0010834, 0.0010864\nuₘ-u, uₙ-u = 0.0010880, 0.0010910\nuₘ-u, uₙ-u = 0.0010923, 0.0010952\nuₘ-u, uₙ-u = 0.0010963, 0.0010992\nuₘ-u, uₙ-u = 0.0011000, 0.0011028\nuₘ-u, uₙ-u = 0.0011034, 0.0011062\nuₘ-u, uₙ-u = 0.0011065, 0.0011093\nuₘ-u, uₙ-u = 0.0011093, 0.0011121\nuₘ-u, uₙ-u = 0.0011119, 0.0011146\nuₘ-u, uₙ-u = 0.0011143, 0.0011169\nuₘ-u, uₙ-u = 0.0011163, 0.0011190\nuₘ-u, uₙ-u = 0.0011182, 0.0011208\nuₘ-u, uₙ-u = 0.0011198, 0.0011224\nuₘ-u, uₙ-u = 0.0011212, 0.0011238\nuₘ-u, uₙ-u = 0.0011224, 0.0011249\nuₘ-u, uₙ-u = 0.0011234, 0.0011259\nuₘ-u, uₙ-u = 0.0011241, 0.0011266\nuₘ-u, uₙ-u = 0.0011247, 0.0011271\nuₘ-u, uₙ-u = 0.0011251, 0.0011275\nuₘ-u, uₙ-u = 0.0011253, 0.0011276\nuₘ-u, uₙ-u = 0.0011253, 0.0011276\nuₘ-u, uₙ-u = 0.0011251, 0.0011274\nuₘ-u, uₙ-u = 0.0011247, 0.0011270\nuₘ-u, uₙ-u = 0.0011242, 0.0011265\nuₘ-u, uₙ-u = 0.0011236, 0.0011258\nuₘ-u, uₙ-u = 0.0011227, 0.0011249\nuₘ-u, uₙ-u = 0.0011217, 0.0011239\nuₘ-u, uₙ-u = 0.0011206, 0.0011228\nuₘ-u, uₙ-u = 0.0011193, 0.0011215\nuₘ-u, uₙ-u = 0.0011179, 0.0011200\nuₘ-u, uₙ-u = 0.0011164, 0.0011184\nuₘ-u, uₙ-u = 0.0011147, 0.0011167\nuₘ-u, uₙ-u = 0.0011128, 0.0011149\nuₘ-u, uₙ-u = 0.0011109, 0.0011129\nuₘ-u, uₙ-u = 0.0011088, 0.0011108\nuₘ-u, uₙ-u = 0.0011067, 0.0011086\nuₘ-u, uₙ-u = 0.0011044, 0.0011063\nuₘ-u, uₙ-u = 0.0011019, 0.0011039\nuₘ-u, uₙ-u = 0.0010994, 0.0011013\nuₘ-u, uₙ-u = 0.0010968, 0.0010987\nuₘ-u, uₙ-u = 0.0010941, 0.0010959\nuₘ-u, uₙ-u = 0.0010912, 0.0010931\nuₘ-u, uₙ-u = 0.0010883, 0.0010901\nuₘ-u, uₙ-u = 0.0010853, 0.0010871\nuₘ-u, uₙ-u = 0.0010822, 0.0010839\nuₘ-u, uₙ-u = 0.0010790, 0.0010807\nuₘ-u, uₙ-u = 0.0010757, 0.0010774\nuₘ-u, uₙ-u = 0.0010723, 0.0010740\nuₘ-u, uₙ-u = 0.0010689, 0.0010706\nuₘ-u, uₙ-u = 0.0010654, 0.0010670\nuₘ-u, uₙ-u = 0.0010618, 0.0010634\nuₘ-u, uₙ-u = 0.0010581, 0.0010597\nuₘ-u, uₙ-u = 0.0010544, 0.0010560\nuₘ-u, uₙ-u = 0.0010506, 0.0010521\nuₘ-u, uₙ-u = 0.0010467, 0.0010482\nuₘ-u, uₙ-u = 0.0010427, 0.0010443\nuₘ-u, uₙ-u = 0.0010387, 0.0010403\nuₘ-u, uₙ-u = 0.0010347, 0.0010362\nuₘ-u, uₙ-u = 0.0010306, 0.0010321\nuₘ-u, uₙ-u = 0.0010264, 0.0010279\nuₘ-u, uₙ-u = 0.0010222, 0.0010236\nuₘ-u, uₙ-u = 0.0010179, 0.0010193\nuₘ-u, uₙ-u = 0.0010136, 0.0010150\nuₘ-u, uₙ-u = 0.0010092, 0.0010106\nuₘ-u, uₙ-u = 0.0010048, 0.0010062\nuₘ-u, uₙ-u = 0.0010003, 0.0010017\nuₘ-u, uₙ-u = 0.0009958, 0.0009972\nuₘ-u, uₙ-u = 0.0009913, 0.0009926\nuₘ-u, uₙ-u = 0.0009867, 0.0009880\nuₘ-u, uₙ-u = 0.0009821, 0.0009834\nuₘ-u, uₙ-u = 0.0009774, 0.0009787\nuₘ-u, uₙ-u = 0.0009727, 0.0009740\nuₘ-u, uₙ-u = 0.0009680, 0.0009693\nuₘ-u, uₙ-u = 0.0009633, 0.0009645\nuₘ-u, uₙ-u = 0.0009585, 0.0009597\nuₘ-u, uₙ-u = 0.0009537, 0.0009549\nuₘ-u, uₙ-u = 0.0009489, 0.0009500\nuₘ-u, uₙ-u = 0.0009440, 0.0009452\nuₘ-u, uₙ-u = 0.0009391, 0.0009403\nuₘ-u, uₙ-u = 0.0009342, 0.0009354\nuₘ-u, uₙ-u = 0.0009293, 0.0009304\nuₘ-u, uₙ-u = 0.0009243, 0.0009255\nuₘ-u, uₙ-u = 0.0009194, 0.0009205\nuₘ-u, uₙ-u = 0.0009144, 0.0009155\nuₘ-u, uₙ-u = 0.0009094, 0.0009105\nuₘ-u, uₙ-u = 0.0009044, 0.0009054\nuₘ-u, uₙ-u = 0.0008993, 0.0009004\nuₘ-u, uₙ-u = 0.0008943, 0.0008953\nuₘ-u, uₙ-u = 0.0008892, 0.0008903\nuₘ-u, uₙ-u = 0.0008842, 0.0008852\nuₘ-u, uₙ-u = 0.0008791, 0.0008801\nuₘ-u, uₙ-u = 0.0008740, 0.0008750\nuₘ-u, uₙ-u = 0.0008689, 0.0008699\nuₘ-u, uₙ-u = 0.0008638, 0.0008648\nuₘ-u, uₙ-u = 0.0008587, 0.0008597\nuₘ-u, uₙ-u = 0.0008536, 0.0008546\nuₘ-u, uₙ-u = 0.0008485, 0.0008495\nuₘ-u, uₙ-u = 0.0008434, 0.0008443\nuₘ-u, uₙ-u = 0.0008383, 0.0008392\nuₘ-u, uₙ-u = 0.0008332, 0.0008341\nuₘ-u, uₙ-u = 0.0008280, 0.0008289\nuₘ-u, uₙ-u = 0.0008229, 0.0008238\nuₘ-u, uₙ-u = 0.0008178, 0.0008187\nuₘ-u, uₙ-u = 0.0008127, 0.0008136\nuₘ-u, uₙ-u = 0.0008076, 0.0008084\nuₘ-u, uₙ-u = 0.0008025, 0.0008033\nuₘ-u, uₙ-u = 0.0007974, 0.0007982\nuₘ-u, uₙ-u = 0.0007923, 0.0007931\nuₘ-u, uₙ-u = 0.0007872, 0.0007880\nuₘ-u, uₙ-u = 0.0007821, 0.0007829\nuₘ-u, uₙ-u = 0.0007770, 0.0007778\nuₘ-u, uₙ-u = 0.0007719, 0.0007727\nuₘ-u, uₙ-u = 0.0007668, 0.0007676\nuₘ-u, uₙ-u = 0.0007618, 0.0007625\nuₘ-u, uₙ-u = 0.0007567, 0.0007575\nuₘ-u, uₙ-u = 0.0007517, 0.0007524\nuₘ-u, uₙ-u = 0.0007466, 0.0007474\nuₘ-u, uₙ-u = 0.0007416, 0.0007423\nuₘ-u, uₙ-u = 0.0007366, 0.0007373\nuₘ-u, uₙ-u = 0.0007316, 0.0007323\nuₘ-u, uₙ-u = 0.0007266, 0.0007273\nuₘ-u, uₙ-u = 0.0007216, 0.0007223\nuₘ-u, uₙ-u = 0.0007167, 0.0007173\nuₘ-u, uₙ-u = 0.0007117, 0.0007124\nuₘ-u, uₙ-u = 0.0007068, 0.0007074\nuₘ-u, uₙ-u = 0.0007018, 0.0007025\nuₘ-u, uₙ-u = 0.0006969, 0.0006976\nuₘ-u, uₙ-u = 0.0006920, 0.0006927\nuₘ-u, uₙ-u = 0.0006871, 0.0006878\nuₘ-u, uₙ-u = 0.0006823, 0.0006829\nuₘ-u, uₙ-u = 0.0006774, 0.0006781\nuₘ-u, uₙ-u = 0.0006726, 0.0006732\nuₘ-u, uₙ-u = 0.0006678, 0.0006684\nuₘ-u, uₙ-u = 0.0006630, 0.0006636\nuₘ-u, uₙ-u = 0.0006582, 0.0006588\nuₘ-u, uₙ-u = 0.0006534, 0.0006540\nuₘ-u, uₙ-u = 0.0006487, 0.0006493\nuₘ-u, uₙ-u = 0.0006440, 0.0006445\nuₘ-u, uₙ-u = 0.0006392, 0.0006398\nuₘ-u, uₙ-u = 0.0006346, 0.0006351\nuₘ-u, uₙ-u = 0.0006299, 0.0006304\nuₘ-u, uₙ-u = 0.0006252, 0.0006258\nuₘ-u, uₙ-u = 0.0006206, 0.0006211\nuₘ-u, uₙ-u = 0.0006160, 0.0006165\nuₘ-u, uₙ-u = 0.0006114, 0.0006119\nuₘ-u, uₙ-u = 0.0006068, 0.0006073\nuₘ-u, uₙ-u = 0.0006022, 0.0006028\nuₘ-u, uₙ-u = 0.0005977, 0.0005982\nuₘ-u, uₙ-u = 0.0005932, 0.0005937\nuₘ-u, uₙ-u = 0.0005887, 0.0005892\nuₘ-u, uₙ-u = 0.0005842, 0.0005847\nuₘ-u, uₙ-u = 0.0005798, 0.0005803\nuₘ-u, uₙ-u = 0.0005753, 0.0005758\nuₘ-u, uₙ-u = 0.0005709, 0.0005714\nuₘ-u, uₙ-u = 0.0005666, 0.0005670\nuₘ-u, uₙ-u = 0.0005622, 0.0005626\nuₘ-u, uₙ-u = 0.0005578, 0.0005583\nuₘ-u, uₙ-u = 0.0005535, 0.0005540\nuₘ-u, uₙ-u = 0.0005492, 0.0005497\nuₘ-u, uₙ-u = 0.0005449, 0.0005454\nuₘ-u, uₙ-u = 0.0005407, 0.0005411\nuₘ-u, uₙ-u = 0.0005364, 0.0005369\nuₘ-u, uₙ-u = 0.0005322, 0.0005327\nuₘ-u, uₙ-u = 0.0005280, 0.0005285\nuₘ-u, uₙ-u = 0.0005239, 0.0005243\nuₘ-u, uₙ-u = 0.0005197, 0.0005201\nuₘ-u, uₙ-u = 0.0005156, 0.0005160\nuₘ-u, uₙ-u = 0.0005115, 0.0005119\nuₘ-u, uₙ-u = 0.0005074, 0.0005078\nuₘ-u, uₙ-u = 0.0005034, 0.0005038\nuₘ-u, uₙ-u = 0.0004994, 0.0004997\nuₘ-u, uₙ-u = 0.0004954, 0.0004957\nuₘ-u, uₙ-u = 0.0004914, 0.0004917\nuₘ-u, uₙ-u = 0.0004874, 0.0004878\nuₘ-u, uₙ-u = 0.0004835, 0.0004838\nuₘ-u, uₙ-u = 0.0004796, 0.0004799\nuₘ-u, uₙ-u = 0.0004757, 0.0004760\nuₘ-u, uₙ-u = 0.0004718, 0.0004722\nuₘ-u, uₙ-u = 0.0004680, 0.0004683\nuₘ-u, uₙ-u = 0.0004641, 0.0004645\nuₘ-u, uₙ-u = 0.0004604, 0.0004607\nuₘ-u, uₙ-u = 0.0004566, 0.0004569\nuₘ-u, uₙ-u = 0.0004528, 0.0004532\nuₘ-u, uₙ-u = 0.0004491, 0.0004494\nuₘ-u, uₙ-u = 0.0004454, 0.0004457\nuₘ-u, uₙ-u = 0.0004417, 0.0004420\nuₘ-u, uₙ-u = 0.0004381, 0.0004384\nuₘ-u, uₙ-u = 0.0004344, 0.0004347\nuₘ-u, uₙ-u = 0.0004308, 0.0004311\nuₘ-u, uₙ-u = 0.0004272, 0.0004275\nuₘ-u, uₙ-u = 0.0004237, 0.0004240\nuₘ-u, uₙ-u = 0.0004201, 0.0004204\nuₘ-u, uₙ-u = 0.0004166, 0.0004169\nuₘ-u, uₙ-u = 0.0004131, 0.0004134\nuₘ-u, uₙ-u = 0.0004096, 0.0004099\nuₘ-u, uₙ-u = 0.0004062, 0.0004065\nuₘ-u, uₙ-u = 0.0004028, 0.0004030\nuₘ-u, uₙ-u = 0.0003993, 0.0003996\nuₘ-u, uₙ-u = 0.0003960, 0.0003962\nuₘ-u, uₙ-u = 0.0003926, 0.0003929\nuₘ-u, uₙ-u = 0.0003893, 0.0003895\nuₘ-u, uₙ-u = 0.0003860, 0.0003862\nuₘ-u, uₙ-u = 0.0003827, 0.0003829\nuₘ-u, uₙ-u = 0.0003794, 0.0003797\nuₘ-u, uₙ-u = 0.0003762, 0.0003764\nuₘ-u, uₙ-u = 0.0003729, 0.0003732\nuₘ-u, uₙ-u = 0.0003697, 0.0003700\nuₘ-u, uₙ-u = 0.0003666, 0.0003668\nuₘ-u, uₙ-u = 0.0003634, 0.0003636\nuₘ-u, uₙ-u = 0.0003603, 0.0003605\nuₘ-u, uₙ-u = 0.0003572, 0.0003574\nuₘ-u, uₙ-u = 0.0003541, 0.0003543\nuₘ-u, uₙ-u = 0.0003510, 0.0003512\nuₘ-u, uₙ-u = 0.0003480, 0.0003482\nuₘ-u, uₙ-u = 0.0003449, 0.0003452\nuₘ-u, uₙ-u = 0.0003419, 0.0003422\nuₘ-u, uₙ-u = 0.0003390, 0.0003392\nuₘ-u, uₙ-u = 0.0003360, 0.0003362\nuₘ-u, uₙ-u = 0.0003331, 0.0003333\nuₘ-u, uₙ-u = 0.0003302, 0.0003304\nuₘ-u, uₙ-u = 0.0003273, 0.0003275\nuₘ-u, uₙ-u = 0.0003244, 0.0003246\nuₘ-u, uₙ-u = 0.0003215, 0.0003217\nuₘ-u, uₙ-u = 0.0003187, 0.0003189\nuₘ-u, uₙ-u = 0.0003159, 0.0003161\nuₘ-u, uₙ-u = 0.0003131, 0.0003133\nuₘ-u, uₙ-u = 0.0003103, 0.0003105\nuₘ-u, uₙ-u = 0.0003076, 0.0003078\nuₘ-u, uₙ-u = 0.0003049, 0.0003051\nuₘ-u, uₙ-u = 0.0003022, 0.0003023\nuₘ-u, uₙ-u = 0.0002995, 0.0002997\nuₘ-u, uₙ-u = 0.0002968, 0.0002970\nuₘ-u, uₙ-u = 0.0002942, 0.0002943\nuₘ-u, uₙ-u = 0.0002915, 0.0002917\nuₘ-u, uₙ-u = 0.0002889, 0.0002891\nuₘ-u, uₙ-u = 0.0002864, 0.0002865\nuₘ-u, uₙ-u = 0.0002838, 0.0002840\nuₘ-u, uₙ-u = 0.0002812, 0.0002814\nuₘ-u, uₙ-u = 0.0002787, 0.0002789\nuₘ-u, uₙ-u = 0.0002762, 0.0002764\nuₘ-u, uₙ-u = 0.0002737, 0.0002739\nuₘ-u, uₙ-u = 0.0002713, 0.0002714\nuₘ-u, uₙ-u = 0.0002688, 0.0002690\nuₘ-u, uₙ-u = 0.0002664, 0.0002665\nuₘ-u, uₙ-u = 0.0002640, 0.0002641\nuₘ-u, uₙ-u = 0.0002616, 0.0002617\nuₘ-u, uₙ-u = 0.0002592, 0.0002594\nuₘ-u, uₙ-u = 0.0002569, 0.0002570\nuₘ-u, uₙ-u = 0.0002545, 0.0002547\nuₘ-u, uₙ-u = 0.0002522, 0.0002524\nuₘ-u, uₙ-u = 0.0002499, 0.0002501\nuₘ-u, uₙ-u = 0.0002477, 0.0002478\nuₘ-u, uₙ-u = 0.0002454, 0.0002455\nuₘ-u, uₙ-u = 0.0002432, 0.0002433\nuₘ-u, uₙ-u = 0.0002409, 0.0002411\nuₘ-u, uₙ-u = 0.0002387, 0.0002389\nuₘ-u, uₙ-u = 0.0002365, 0.0002367\nuₘ-u, uₙ-u = 0.0002344, 0.0002345\nuₘ-u, uₙ-u = 0.0002322, 0.0002323\nuₘ-u, uₙ-u = 0.0002301, 0.0002302\nuₘ-u, uₙ-u = 0.0002280, 0.0002281\nuₘ-u, uₙ-u = 0.0002259, 0.0002260\nuₘ-u, uₙ-u = 0.0002238, 0.0002239\nuₘ-u, uₙ-u = 0.0002217, 0.0002218\nuₘ-u, uₙ-u = 0.0002197, 0.0002198\nuₘ-u, uₙ-u = 0.0002176, 0.0002178\nuₘ-u, uₙ-u = 0.0002156, 0.0002157\nuₘ-u, uₙ-u = 0.0002136, 0.0002137\nuₘ-u, uₙ-u = 0.0002116, 0.0002118\nuₘ-u, uₙ-u = 0.0002097, 0.0002098\nuₘ-u, uₙ-u = 0.0002077, 0.0002078\nuₘ-u, uₙ-u = 0.0002058, 0.0002059\nuₘ-u, uₙ-u = 0.0002039, 0.0002040\nuₘ-u, uₙ-u = 0.0002020, 0.0002021\nuₘ-u, uₙ-u = 0.0002001, 0.0002002\nuₘ-u, uₙ-u = 0.0001982, 0.0001983\nuₘ-u, uₙ-u = 0.0001964, 0.0001965\nuₘ-u, uₙ-u = 0.0001945, 0.0001946\nuₘ-u, uₙ-u = 0.0001927, 0.0001928\nuₘ-u, uₙ-u = 0.0001909, 0.0001910\nuₘ-u, uₙ-u = 0.0001891, 0.0001892\nuₘ-u, uₙ-u = 0.0001873, 0.0001874\nuₘ-u, uₙ-u = 0.0001856, 0.0001857\nuₘ-u, uₙ-u = 0.0001838, 0.0001839\nuₘ-u, uₙ-u = 0.0001821, 0.0001822\nuₘ-u, uₙ-u = 0.0001804, 0.0001805\nuₘ-u, uₙ-u = 0.0001787, 0.0001788\nuₘ-u, uₙ-u = 0.0001770, 0.0001771\nuₘ-u, uₙ-u = 0.0001753, 0.0001754\nuₘ-u, uₙ-u = 0.0001737, 0.0001738\nuₘ-u, uₙ-u = 0.0001720, 0.0001721\nuₘ-u, uₙ-u = 0.0001704, 0.0001705\nuₘ-u, uₙ-u = 0.0001688, 0.0001689\nuₘ-u, uₙ-u = 0.0001672, 0.0001673\nuₘ-u, uₙ-u = 0.0001656, 0.0001657\nuₘ-u, uₙ-u = 0.0001640, 0.0001641\nuₘ-u, uₙ-u = 0.0001624, 0.0001625\nuₘ-u, uₙ-u = 0.0001609, 0.0001610\nuₘ-u, uₙ-u = 0.0001594, 0.0001594\nuₘ-u, uₙ-u = 0.0001579, 0.0001579\nuₘ-u, uₙ-u = 0.0001563, 0.0001564\nuₘ-u, uₙ-u = 0.0001549, 0.0001549\nuₘ-u, uₙ-u = 0.0001534, 0.0001534\nuₘ-u, uₙ-u = 0.0001519, 0.0001520\nuₘ-u, uₙ-u = 0.0001504, 0.0001505\nuₘ-u, uₙ-u = 0.0001490, 0.0001491\nuₘ-u, uₙ-u = 0.0001476, 0.0001477\nuₘ-u, uₙ-u = 0.0001462, 0.0001462\nuₘ-u, uₙ-u = 0.0001448, 0.0001448\nuₘ-u, uₙ-u = 0.0001434, 0.0001434\nuₘ-u, uₙ-u = 0.0001420, 0.0001421\nuₘ-u, uₙ-u = 0.0001406, 0.0001407\nuₘ-u, uₙ-u = 0.0001393, 0.0001393\nuₘ-u, uₙ-u = 0.0001379, 0.0001380\nuₘ-u, uₙ-u = 0.0001366, 0.0001367\nuₘ-u, uₙ-u = 0.0001353, 0.0001353\nuₘ-u, uₙ-u = 0.0001340, 0.0001340\nuₘ-u, uₙ-u = 0.0001327, 0.0001327\nuₘ-u, uₙ-u = 0.0001314, 0.0001315\nuₘ-u, uₙ-u = 0.0001301, 0.0001302\nuₘ-u, uₙ-u = 0.0001289, 0.0001289\nuₘ-u, uₙ-u = 0.0001276, 0.0001277\nuₘ-u, uₙ-u = 0.0001264, 0.0001264\nuₘ-u, uₙ-u = 0.0001251, 0.0001252\nuₘ-u, uₙ-u = 0.0001239, 0.0001240\nuₘ-u, uₙ-u = 0.0001227, 0.0001228\nuₘ-u, uₙ-u = 0.0001215, 0.0001216\nuₘ-u, uₙ-u = 0.0001203, 0.0001204\nuₘ-u, uₙ-u = 0.0001192, 0.0001192\nuₘ-u, uₙ-u = 0.0001180, 0.0001181\nuₘ-u, uₙ-u = 0.0001169, 0.0001169\nuₘ-u, uₙ-u = 0.0001157, 0.0001158\nuₘ-u, uₙ-u = 0.0001146, 0.0001146\nuₘ-u, uₙ-u = 0.0001135, 0.0001135\nuₘ-u, uₙ-u = 0.0001124, 0.0001124\nuₘ-u, uₙ-u = 0.0001113, 0.0001113\nuₘ-u, uₙ-u = 0.0001102, 0.0001102\nuₘ-u, uₙ-u = 0.0001091, 0.0001091\nuₘ-u, uₙ-u = 0.0001080, 0.0001081\n"
!     }
! ]
