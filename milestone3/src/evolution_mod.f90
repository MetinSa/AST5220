 module evolution_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use rec_mod
  implicit none

  ! Accuracy parameters
  real(dp),     parameter, private :: a_init   = 1.d-8
  real(dp),     parameter, private :: x_init   = log(a_init)
  real(dp),     parameter, private :: k_min    = 0.1d0 * H_0 / c
  real(dp),     parameter, private :: k_max    = 1.d3  * H_0 / c
  integer(i4b), parameter          :: n_k      = 100
  integer(i4b), parameter, private :: lmax_int = 6

  ! Perturbation quantities
  real(dp), allocatable, dimension(:,:,:) :: Theta
  real(dp), allocatable, dimension(:,:)   :: delta
  real(dp), allocatable, dimension(:,:)   :: delta_b
  real(dp), allocatable, dimension(:,:)   :: Phi
  real(dp), allocatable, dimension(:,:)   :: Psi
  real(dp), allocatable, dimension(:,:)   :: v
  real(dp), allocatable, dimension(:,:)   :: v_b
  real(dp), allocatable, dimension(:,:)   :: dPhi
  real(dp), allocatable, dimension(:,:)   :: dPsi
  real(dp), allocatable, dimension(:,:)   :: dv_b
  real(dp), allocatable, dimension(:,:,:) :: dTheta

  ! Fourier mode list
  real(dp), allocatable, dimension(:) :: ks

  ! Book-keeping variables
  real(dp),     private :: k_current
  real(dp), allocatable, dimension(:) :: derivative
  real(dp), allocatable, dimension(:) :: k_plot             ! k-values which we plot
  integer(i4b), private :: npar = 6+lmax_int

contains


  ! NB!!! New routine for 4th milestone only; disregard until then!!!
  subroutine get_hires_source_function(k, x, S)
    implicit none

    real(dp), pointer, dimension(:),   intent(out) :: k, x
    real(dp), pointer, dimension(:,:), intent(out) :: S

    integer(i4b) :: i, j
    real(dp)     :: g, dg, ddg, tau, dt, ddt, H_p, dH_p, ddHH_p, Pi, dPi, ddPi
    real(dp), allocatable, dimension(:,:) :: S_lores

    ! Task: Output a pre-computed 2D array (over k and x) for the 
    !       source function, S(k,x). Remember to set up (and allocate) output 
    !       k and x arrays too. 	
    !
    ! Substeps:
    !   1) First compute the source function over the existing k and x
    !      grids
    !   2) Then spline this function with a 2D spline
    !   3) Finally, resample the source function on a high-resolution uniform
    !      5000 x 5000 grid and return this, together with corresponding
    !      high-resolution k and x arrays

  end subroutine get_hires_source_function



  !-----------------------------------------------------------------
  ! Milestone 3 subroutines
  !-----------------------------------------------------------------

  ! Routine for initializing and solving the Boltzmann and Einstein equations
  subroutine initialize_perturbation_eqns
    implicit none

    integer(i4b) :: l, i
    real(dp)     :: cksH_p, dtau

    ! Allocate arrays for perturbation quantities
    allocate(Theta(0:n_t, 0:lmax_int, n_k))
    allocate(delta(0:n_t, n_k))
    allocate(delta_b(0:n_t, n_k))
    allocate(v(0:n_t, n_k))
    allocate(v_b(0:n_t, n_k))
    allocate(Phi(0:n_t, n_k))
    allocate(Psi(0:n_t, n_k))
    allocate(dPhi(0:n_t, n_k))
    allocate(dPsi(0:n_t, n_k))
    allocate(dv_b(0:n_t, n_k))
    allocate(dTheta(0:n_t, 0:lmax_int, n_k))
    allocate(ks(n_k))

    ! Task: Initialize k-grid, ks; quadratic between k_min and k_max
    do i = 1, n_k
       ks(i) = k_min + (k_max - k_min)*((i-1.d0)/(n_k-1.d0))**2
    end do

    ! Task: Set up initial conditions for the Boltzmann and Einstein equations
    Phi(0,:)     = 1.d0
    delta(0,:)   = 1.5d0 * Phi(0,:)
    delta_b(0,:) = delta(0,:)
    Theta(0,0,:) = 0.5d0 * Phi(0,:)
       
    do i = 1, n_k
       ! Substitution (efficiency)
       cksH_p = c*ks(i)/get_H_p(x_init)
       dtau = get_dtau(x_init)

       v(0,i)       = (cksH_p*Phi(0,i))/2.d0
       v_b(0,i)     = v(0,i)
       Theta(0,1,i) = -(cksH_p*Phi(0,i))/6.d0
       Theta(0,2,i) =  -(20.d0*cksH_p*Theta(0,1,i))/(45.d0*dtau) 
       do l = 3, lmax_int
          Theta(0,l,i) = - (l/(2.d0*l + 1.d0)) * ((cksH_p*Theta(0,l-1,i))/dtau)
       end do

    end do

  end subroutine initialize_perturbation_eqns


  ! Routine which integrates and solves the Boltzmann and Einstein equations
  subroutine integrate_perturbation_eqns
    implicit none

    integer(i4b) :: i, i_tc, j, k, l
    real(dp)     :: x_init, x_tc
    real(dp)     :: eps, hmin, h1, x, a, H_p, ckH_p, dtau, ck

    real(dp), allocatable, dimension(:) :: y, y_tight_coupling

    eps    = 1.d-8
    hmin   = 0.d0
    h1     = 1.d-5

    allocate(y(npar))
    allocate(derivative(npar))
    allocate(y_tight_coupling(7))
    allocate(k_plot(6))

    ! The different k's in the interval of intereset which we will plot
    k_plot(1)   = 1
    k_plot(2)   = 10
    k_plot(3)   = 30
    k_plot(4)   = 50
    k_plot(5)   = 80
    k_plot(6)   = 100

    ! Opening files which we will write to
    open(1, file = 'milestone3.dat', action = 'write', status = 'replace') 

    ! Propagate each k-mode independently
    do k = 1, n_k

       k_current = ks(k)  ! Store k_current as a global module variable
       ck = c*k_current   ! Substitution to increase efficiency

       ! Initialize equation set for tight coupling
       y_tight_coupling(1) = delta(0,k)
       y_tight_coupling(2) = delta_b(0,k)
       y_tight_coupling(3) = v(0,k)
       y_tight_coupling(4) = v_b(0,k)
       y_tight_coupling(5) = Phi(0,k)
       y_tight_coupling(6) = Theta(0,0,k)
       y_tight_coupling(7) = Theta(0,1,k)
       
       ! Find the time to which tight coupling is assumed, 
       ! and integrate equations to that time
       x_tc = get_tight_coupling_time(k_current)

       ! Task: Integrate from x_init until the end of tight coupling, using
       !       the tight coupling equations

       ! Iteration in tight coupled regime
       i_tc = 1

       ! x_t is from x_init til x_today, defined in time_mod
       do while (x_t(i_tc+1) < x_tc)

          !-----------------------------------
          ! Pre tight coupling
          !----------------------------------

          ! Substitutions
          x = x_t(i_tc+1)
          a = exp(x)
          H_p = get_H_p(x)
          ckH_p = ck/H_p
          dtau = get_dtau(x)

          ! Compute derivatives
          call odeint(y_tight_coupling, x_t(i_tc), x, eps, h1, hmin, dy_tcdx, bsstep, output)

          ! Saving the integrated values
          delta(i_tc,k)    = y_tight_coupling(1)
          delta_b(i_tc,k)  = y_tight_coupling(2)
          v(i_tc,k)        = y_tight_coupling(3)
          v_b(i_tc,k)      = y_tight_coupling(4)
          Phi(i_tc,k)      = y_tight_coupling(5)
          Theta(i_tc,0,k)  = y_tight_coupling(6)
          Theta(i_tc,1,k)  = y_tight_coupling(7)
          Theta(i_tc,2,k)  = -((20.d0*ckH_p)/(45.d0*dtau))*Theta(i_tc,1,k)
          Psi(i_tc,k)     = -Phi(i_tc,k) - 12.d0*Omega_r*Theta(i_tc,2,k)*(H_0/(ck*a))**2


          do l = 3, lmax_int
             Theta(i_tc,l,k) = -(l/(2.d0*l + 1.d0))*(ckH_p/dtau)*Theta(i_tc,l-1,k)
          end do

          ! Save derivatives
          call dy_tcdx(x, y_tight_coupling, derivative)
          dv_b(i_tc,k)       = derivative(4)
          dPhi(i_tc,k)       = derivative(5)
          dTheta(i_tc,0,k)   = derivative(6)
          dTheta(i_tc,1,k)   = derivative(7)
          dTheta(i_tc,2,k)   = (2.d0/5.d0)*ckH_p*Theta(i_tc,1,k) - (3.d0/5.d0)*ckH_p*Theta(i_tc,3,k) + dtau*0.9d0*Theta(i_tc,2,k)   
          dPsi(i_tc,k)     = -dPhi(i_tc,k) - ((12.d0*H_0**2)/((ck*a)**2))*Omega_r*(-2.d0*Theta(i_tc,2,k) + dTheta(i_tc,2,k))

          do l = 3, lmax_int-1
             dTheta(i_tc,l,k) = (l*ckH_p*Theta(i_tc,l-1,k))/(2.d0*l + 1.d0) - ((l+1.d0)*ckH_p*Theta(i_tc,l+1,k))/(2.d0*l + 1.d0) + dtau*Theta(i_tc,l,k)
          end do

          ! Writing to file for pre tc
          do j = 1, 6
             if (k == k_plot(j)) then
                write(1, '(11(E17.8))') x_t(i_tc), delta(i_tc,k), delta_b(i_tc,k), v(i_tc,k), v_b(i_tc,k), Phi(i_tc,k), Psi(i_tc,k), dPhi(i_tc,k), dPsi(i_tc,k), Theta(i_tc,0,k), Theta(i_tc,1,k)
             end if
          end do


          ! Increasing iteration
          i_tc = i_tc +1
          

       end do

       !-----------------------------------
       ! Post tight coupling
       !----------------------------------

       ! Task: Set up variables for integration from the end of tight coupling 
       ! until today

       ! Initializing set of equations y, with the final y_tight_coupling values
       y(1:7) = y_tight_coupling(1:7)
       y(8)   = Theta(i_tc-1,2,k)

       do l = 3, lmax_int
          y(6+l) = Theta(i_tc-1,l,k)
       end do

       do i = i_tc, n_t
          ! Task: Integrate equations from tight coupling to today

          ! Substitutions
          x = x_t(i)
          a = exp(x)

          ! Compute derivatives
          call odeint(y, x_t(i-1), x, eps, h1, hmin, dydx, bsstep, output)

          ! Task: Store variables at time step i in global variables
          delta(i,k)   = y(1)
          delta_b(i,k) = y(2)
          v(i,k)       = y(3)
          v_b(i,k)     = y(4)
          Phi(i,k)     = y(5)

          do l = 0, lmax_int
             Theta(i,l,k) = y(6+l)
          end do

          Psi(i,k)     = -Phi(i,k) - 12.d0*Omega_r*Theta(i,2,k)*(H_0/(ck*a))**2

          ! Task: Store derivatives that are required for C_l estimation
          call dydx(x, y, derivative)

          dPhi(i,k)     = derivative(4)
          dv_b(i,k)     = derivative(5)

          do l = 0, lmax_int
             dTheta(i,l,k) = derivative(6+l) 
          end do

          dPsi(i,k)     = -dPhi(i,k) - ((12.d0*H_0**2)/((ck*a)**2))*Omega_r*(-2.d0*Theta(i,2,k) + dTheta(i,2,k))          
          
          do j = 1, 6
             if (k == k_plot(j)) then
                write(1, '(11(E17.8))') x_t(i), delta(i,k), delta_b(i,k), v(i,k), v_b(i,k), Phi(i,k), Psi(i,k), dPhi(i,k), dPsi(i,k), Theta(i,0,k), Theta(i,1,k)
             end if
          end do

       end do

    end do

    deallocate(y_tight_coupling)
    deallocate(y)
    deallocate(derivative)


  end subroutine integrate_perturbation_eqns


  ! Computing the derivatives of y_tight_coupling set of equations
  subroutine dy_tcdx(x, y_tc, derivative)
    implicit none

    real(dp),                intent(in)  :: x
    real(dp),  dimension(:), intent(in)  :: y_tc
    real(dp),  dimension(:), intent(out) :: derivative

    real(dp)  :: delta, delta_b, v, v_b, Phi, Psi, Theta0, Theta1, Theta2
    real(dp)  :: ddelta, ddelta_b, dv, dv_b, dPhi, dTheta0, dTheta1, dTheta2
    real(dp)  :: R, q, a, H_p, dH_p, dtau, ddtau, ckH_p

    ! Reading in initial values
    delta         = y_tc(1)
    delta_b       = y_tc(2)
    v             = y_tc(3)
    v_b           = y_tc(4)
    Phi           = y_tc(5)
    Theta0        = y_tc(6)
    Theta1        = y_tc(7)

    ! Defining some quantities and substitutions
    a             = exp(x)
    H_p           = get_H_p(x)
    dH_p          = get_dH_p(x)
    dtau          = get_dtau(x)
    ddtau         = get_ddtau(x)
    ckH_p         = c*k_current/H_p


    ! Computing the derivatives and the quantities contained in their expressions
    R             = (4.d0*Omega_r)/(3.d0*Omega_b*a)
    Theta2        = -((20.d0*ckH_p)/(45.d0*dtau))*Theta1
    Psi           = -Phi - 12.d0*Omega_r*Theta2*(H_0/(c*k_current*a))**2
    dPhi          = Psi - ((ckH_p**2)/3.d0)*Phi + (((H_0/H_p)**2)/2.d0)*( ((Omega_m*delta)/a) + ((Omega_b*delta_b)/a) + ((4.d0*Omega_r*Theta0)/a**2) )
    ddelta        = ckH_p*v - 3.d0*dPhi
    ddelta_b      = ckH_p*v_b - 3.d0*dPhi
    dv            = -v - ckH_p*Psi
    dTheta0       = -ckH_p*Theta1 - dPhi
    q             = ( -((1.d0 - 2.d0*R)*dtau + (1.d0 + R)*ddtau)*(3.d0*Theta1 + v_b) - ckH_p*Psi + ((1.d0 - (dH_p/H_p))*ckH_p*(-Theta0 + 2.d0*Theta2)) - ckH_p*dTheta0 )/( (1.d0 + R)*dtau + dH_p/H_p - 1.d0 )
    dv_b          = ( -v_b - ckH_p*Psi + R*(q + ckH_p*(-Theta0 + 2.d0*Theta2) - ckH_p*Psi ))/(1.d0 + R)
    dTheta1       = (q - dv_b)/3.d0

    ! Saving the derivatives 
    derivative(1) = ddelta
    derivative(2) = ddelta_b
    derivative(3) = dv
    derivative(4) = dv_b
    derivative(5) = dPhi
    derivative(6) = dTheta0
    derivative(7) = dTheta1

  end subroutine dy_tcdx



  ! Computing the derivatives of y set of equations
  subroutine dydx(x, y, derivative)
    implicit none

    real(dp),                intent(in)  :: x
    real(dp),  dimension(:), intent(in)  :: y
    real(dp),  dimension(:), intent(out) :: derivative

    integer(i4b) :: l
    real(dp)  :: delta, delta_b, v, v_b, Phi, Psi
    real(dp)  :: Theta0, Theta1, Theta2, Theta3, Theta4, Theta5, Theta6
    real(dp)  :: ddelta, ddelta_b, dv, dv_b, dPhi, dTheta0, dTheta1, dTheta2
    real(dp)  :: R, q, a, H_p, dH_p, dtau, ddtau, ckH_p

    ! Reading in initial values
    delta         = y(1)
    delta_b       = y(2)
    v             = y(3)
    v_b           = y(4)
    Phi           = y(5)
    Theta0        = y(6)
    Theta1        = y(7) 
    Theta2        = y(8) 
    Theta3        = y(9) 
    Theta4        = y(10) 
    Theta5        = y(11) 
    Theta6        = y(12) 
   
    ! Defining some quantities and substitutions
    a             = exp(x)
    H_p           = get_H_p(x)
    dH_p          = get_dH_p(x)
    dtau          = get_dtau(x)
    ddtau         = get_ddtau(x)
    ckH_p         = (c*k_current)/H_p

    ! Computing the derivatives and the quantities contained in their expressions
    R             = (4.d0*Omega_r)/(3.d0*Omega_b*a)
   ! Theta2        = -((20.d0*ckH_p)/(45.d0*dtau))*Theta1
    Psi           = -Phi - 12.d0*Omega_r*Theta2*(H_0/(c*k_current*a))**2
    dPhi          = Psi - ((ckH_p**2)/3.d0)*Phi + (((H_0/H_p)**2)/2.d0)*( ((Omega_m*delta)/a) + ((Omega_b*delta_b)/a) + ((4.d0*Omega_r*Theta0)/a**2) )
    ddelta        = ckH_p*v - 3.d0*dPhi
    ddelta_b      = ckH_p*v_b - 3.d0*dPhi
    dv            = -v - ckH_p*Psi

    ! Manully computing all dTheta but 9,10,11
    dTheta0       = -ckH_p*Theta1 - dPhi
    dv_b          = -v_b - ckH_p*Psi + dtau*R*(3.d0*Theta1 + v_b)
    dTheta1       = ((ckH_p*Theta0)/3.d0) - ((2.d0*ckH_p*Theta2)/3.d0) + ((ckH_p*Psi)/3.d0) + dtau*(Theta1 + (v_b/3.d0))
    dTheta2   = (2.d0/5.d0)*ckH_p*Theta1 - (3.d0/5.d0)*ckH_p*Theta3 + dtau*0.9d0*Theta2    
    do l = 9, 11
       derivative(l) = (l*ckH_p*y(l-1))/(2.d0*l + 1.d0) - ((l+1.d0)*ckH_p*y(l+1))/(2.d0*l + 1.d0) + dtau*y(l)
    end do
    
    derivative(12) = ckH_p*y(11.d0) - (c*(lmax_int +1.d0)*y(12))/(H_p*get_eta(x)) + dtau*y(12)



    ! Saving the derivatives 
    derivative(1) = ddelta
    derivative(2) = ddelta_b
    derivative(3) = dv
    derivative(4) = dv_b
    derivative(5) = dPhi
    derivative(6) = dTheta0
    derivative(7) = dTheta1
    derivative(8) = dTheta2

  end subroutine dydx



  !-----------------------------------------------------------------
  ! Functions
  !-----------------------------------------------------------------

  ! Task: Complete the following routine, such that it returns the time at which
  !       tight coupling ends. In this project, we define this as either when
  !       dtau < 10 or c*k/(H_p*dt) > 0.1 or x > x(start of recombination)

  ! Extracting the time at which tight coupling ends
  function get_tight_coupling_time(k)
    implicit none

    real(dp), intent(in)  :: k
    real(dp)              :: x, get_tight_coupling_time
    integer(i4b)          :: i, n

    ! Setting up semi x-grid from x_init until x_start_rec 
    x = x_init
    ! Checking for 1000 different x values
    n = 1000

    do i = 0, n
       x = x + (x_start_rec - x_init) / (n - 1.d0)
       
       ! Check if TC conditions are fulfilled and extract TC time
       if ((x > x_start_rec) .or. (abs((c*k)/(get_H_p(x)*get_dtau(x))) > 0.1) .or. (abs(get_dtau(x)) < 10.d0)) then

          get_tight_coupling_time = x
          exit

       end if
    end do

  end function get_tight_coupling_time

end module evolution_mod
