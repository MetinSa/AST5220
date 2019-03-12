module rec_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use spline_1D_mod
  implicit none

!-----------------------------------------------------------------------
! Initializing variables and constants
!-----------------------------------------------------------------------


  integer(i4b),                        private :: n                         ! Number of grid points
  real(dp), allocatable, dimension(:), private :: x_rec, a_rec, z_rec       ! x-grid for recombination
  real(dp), allocatable, dimension(:), private :: tau, tau2, tau22          ! tau, tau' and tau''
  real(dp), allocatable, dimension(:), private :: n_e, n_e2                 ! Electron density
  real(dp), allocatable, dimension(:), private :: logn_e, logn_e2           ! (log of) electron density, n_e
  real(dp), allocatable, dimension(:), private :: logtau, logtau2, logtau22 ! (log of) optical depth, tau
  real(dp), allocatable, dimension(:), private :: s_tau, s_tau2, s_tau22    ! Splined quantities of tau
  real(dp), allocatable, dimension(:), private :: g, g2, g22                ! Visibility function g, g' and g''
  real(dp), allocatable, dimension(:), private :: s_g, s_g2, s_g22          ! Splined visibility function

contains

!-----------------------------------------------------------------------
! Initializing the module
!-----------------------------------------------------------------------


  subroutine initialize_rec_mod
    implicit none
    
    integer(i4b) :: i, j, k
    real(dp)     :: saha_limit, y, T_b, n_b, dydx, xmin, xmax, dx_rec, f, n_e0, X_e0, xstart, xstop, gamma, stepsize, stepsize_min, yp1, yp2, eps
    logical(lgt) :: use_saha
    real(dp), allocatable, dimension(:) :: X_e           ! Fractional electron density, n_e / n_H

    saha_limit   = 0.99d0                                ! Switch from Saha to Peebles when X_e < 0.99
    xstart       = log(1.d-10)                           ! Start grids at a = 10^-10
    xstop        = 0.d0                                  ! Stop  grids at a = 1
    n            = 1000                                  ! Number of grid points between xstart and xstopo

    ! Variables and constants connected to ODE-Solver and Spline
    eps          = 1.d-8
    stepsize_min = 1.d-5
    yp1          = 1.d30
    yp2          = 1.d30

    ! Allocating arrays
    allocate(x_rec(n))
    allocate(a_rec(n))
    allocate(z_rec(n))
    allocate(X_e(n))
    allocate(tau(n))
    allocate(tau2(n))
    allocate(tau22(n))
    allocate(n_e(n))
    allocate(n_e2(n))
    allocate(logn_e(n))
    allocate(logn_e2(n))
    allocate(logtau(n))
    allocate(logtau2(n))
    allocate(logtau22(n))
    allocate(g(n))
    allocate(g2(n))
    allocate(g22(n))
    allocate(s_tau(n))
    allocate(s_tau2(n))
    allocate(s_tau22(n))
    allocate(s_g(n))
    allocate(s_g2(n))
    allocate(s_g22(n))

    
    ! Initializing and filling the x_rec array
    x_rec(1) = xstart
    x_rec(n) = xstop
    dx_rec = (xstop - xstart)/(n-1.d0)
    
    do i = 1, n-1
       x_rec(i+1) = x_rec(i) + dx_rec
    end do
    
    ! Computing the corresponding a_rec and z_rec grids for x_rec
    a_rec = exp(x_rec)
    z_rec = 1.d0/a_rec - 1.d0

    ! Stepsize to be used in odeint
    stepsize = abs((x_rec(1) - x_rec(2)))

    ! Starting out by using the Saha equation
    use_saha = .true.
    do i = 1, n

       ! Computing the number density of Baryons
       n_b = (Omega_b*rho_c)/(m_H*a_rec(i)**3.d0)
       
       ! Computing Saha if statement is true
       if (use_saha) then

          ! Substitution and baryon temperature
          T_b = T_0/a_rec(i)
          gamma = (1.d0/(n_b))*((m_e*k_b*T_b)/(2.d0*PI*hbar*hbar))**(1.5d0)*exp(-epsilon_0/(k_b*T_b))
    
          ! Computing the Saha equation ( Electron fraction)
          !X_e(i) =(-gamma + sqrt(gamma**2 + 4.d0*gamma))/2.d0          ! Original quadratic solution
          X_e(i) = 2.d0/(sqrt(1.d0 + 4.d0/gamma)+1.d0)                  ! More stable quadratic solution
  
          ! Computing Peebles equation if statement is false
          if (X_e(i) < saha_limit) use_saha = .false.

       else

          ! Using the odeint to solve Peebles equation along with the dX_edx subroutine
          X_e(i) = X_e(i-1)
          call odeint(X_e(i:i), x_rec(i-1), x_rec(i), eps, stepsize, stepsize_min, dX_edx, rkqs, output) 

       end if

       ! Finaly computing the electron number density
       n_e(i) = X_e(i)*n_b
       
    end do
    

    ! loging and splining the electron density n_e
    logn_e = log(n_e)
    call spline(x_rec, logn_e, yp1, yp2, logn_e2)
    
    ! Reverse integration (initializing)
    tau(n) = 0.d0
    
    do i = n-1 , 1, -1
       
       ! Using odeint to compute the optical depth tau with the use of the dtaudx subroutine
       tau(i) = tau(i+1)
       call odeint(tau(i:i), x_rec(i+1), x_rec(i), eps, stepsize, stepsize_min, dtaudx, rkqs, output)
    
    end do


    tau(n) = tau(n-1)
    ! Computing the splined (log of) optical depth
    ! Computing the splined second derivative of (log of) optical depth
    logtau = log(tau)

    call spline(x_rec, tau, yp1, yp2, tau2)
    call spline(x_rec, tau2, yp1, yp2, tau22)

    call spline(x_rec, logtau, yp1, yp2, logtau2)
    call spline(x_rec, logtau2, yp1, yp2, logtau22)

!-----------------------------------------------------------------------
! Computing the optical depth and visibility function
!-----------------------------------------------------------------------


    ! Computing splined tau, tau' and tau''
    do i = 1, n

       s_tau(i) = get_tau(x_rec(i))
       s_tau2(i) = get_dtau(x_rec(i))
       s_tau22(i) = get_ddtau(x_rec(i))

    end do

    ! Computing the visibility function
    do i = 1, n

       g(i) = -s_tau2(i)*exp(-s_tau(i))

    end do

    ! Computing splined visibility function
    ! Computing splined second derivative of visibility function
    call spline(x_rec, g, yp1, yp2, g2)
    call spline(x_rec, g2, yp1, yp2, g22)

    ! Computing splined g, g' and g''
    do i = 1, n
       
       s_g(i) = get_g(x_rec(i))
       s_g2(i) = get_dg(x_rec(i))
       s_g22(i) = get_ddg(x_rec(i))

    end do
    write(*,*) s_tau22
!-----------------------------------------------------------------------
! Writing data to files
!-----------------------------------------------------------------------


    open(1, file = 'milestone2.dat', action = 'write', status = 'replace')
    
    do i = 1, n
       write(1, "(11(E17.8E3))") x_rec(i), a_rec(i), z_rec(i) , X_e(i), get_n_e(x_rec(i)), s_tau(i), s_tau2(i), s_tau22(i), s_g(i), s_g2(i), s_g22(i)
    end do

    close(1)

  end subroutine initialize_rec_mod

!-----------------------------------------------------------------------
! Subroutines
!-----------------------------------------------------------------------


   subroutine dX_edx(x, X_e, derivative)
     ! Subroutine which computes Peebles equation

     implicit none
     
     real(dp),               intent(in)  :: x
     real(dp), dimension(:), intent(in)  :: X_e
     real(dp), dimension(:), intent(out) :: derivative
     real(dp)                            :: T_b, n_b, phi2, alpha2, beta, beta2, n1s, Lambda_alpha, Lambda, C_r, a, H
     
     ! Hubble parameter and the scale factor
     H            = get_H(x)
     a            = exp(x)

     ! Baryon temperature and number density
     T_b          = T_0/a
     n_b          = (Omega_b*rho_c)/(m_H*a**3)

     ! Substitutions to enter the Peebles Equation
     phi2         = 0.448d0 * log(epsilon_0/(k_b*T_b))
     alpha2       = ((64.d0*PI)/(sqrt(27.d0*PI)))*((alpha/m_e)**2.d0)*sqrt(epsilon_0/(k_b*T_b))*phi2*(hbar*hbar/c)
     beta         = alpha2*((m_e*k_b*T_b)/(2.d0*PI*hbar*hbar))**(1.5d0)*exp(-epsilon_0/(k_b*T_b))
     ! beta2        = beta*exp(3.d0*epsilon_0/(4.d0*k_b*T_b))  ! Original expresssion (less stable?)
     beta2 = alpha2*(m_e*k_b*T_b/(2d0*PI*hbar*hbar))**1.5d0*exp(-1.d0*epsilon_0/(4.d0*k_b*T_b))
     n1s          = (1.d0-X_e(1))*n_b
     Lambda_alpha = H*((3.d0*epsilon_0)**3.d0/(n1s*(8.d0*PI)**2.d0)) / (c*hbar)**3.d0
     Lambda       = 8.227d0
     C_r          = (Lambda + Lambda_alpha)/(Lambda + Lambda_alpha + beta2)

     ! Peebles equation
     derivative = (C_r/H)*(beta*(1.d0-X_e(1)) - n_b*alpha2*X_e(1)**2.d0)

   end subroutine dX_edx


   subroutine dtaudx(x, tau, derivative)
     ! Subroutine which computes the optical depth

     implicit none

     real(dp),                intent(in)  :: x
     real(dp), dimension(:),  intent(in)  :: tau
     real(dp), dimension(:),  intent(out) :: derivative

     ! The optical depth
     derivative = -(get_n_e(x)*sigma_T*c)/(get_H(x))

   end subroutine dtaudx


!-----------------------------------------------------------------------
! Functions
!-----------------------------------------------------------------------


  ! Computing n_e at arbitrary x, using precomputed information
  function get_n_e(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_n_e
    
    ! Exponentiating because we have log n_e
    get_n_e = exp(splint(x_rec, logn_e, logn_e2, x))

  end function get_n_e


  ! Computing tau at arbitrary x, using precomputed information
  function get_tau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_tau

    get_tau = exp(splint(x_rec, logtau, logtau2, x))

  end function get_tau


  ! Computing the derivative of tau at arbitrary x, using precomputed information
  function get_dtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dtau

    get_dtau = get_tau(x)*splint_deriv(x_rec, logtau, logtau2, x)

  end function get_dtau


  ! Computing the second derivative of tau at arbitrary x,  using precomputed information
  function get_ddtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddtau

    get_ddtau = splint(x_rec, logtau2, logtau22, x)*get_tau(x) + get_dtau(x)**2.d0/get_tau(x)

  end function get_ddtau


  ! Computing the visibility function, g, at arbitray x
  function get_g(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_g

    get_g = splint(x_rec, g, g2, x)

  end function get_g


  ! Computing the derivative of the visibility function, g, at arbitray x
  function get_dg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dg

    get_dg = splint_deriv(x_rec, g, g2, x)

  end function get_dg


  ! Computing the second derivative of the visibility function, g, at arbitray x
  function get_ddg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddg

   get_ddg =  splint(x_rec, g2, g22, x)

  end function get_ddg


end module rec_mod
