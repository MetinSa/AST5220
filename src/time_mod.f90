module time_mod
  use healpix_types
  use params
  use spline_1D_mod
  use ode_solver
  implicit none


  !-----------------------------------------------------------------------
  ! Initializing variables and constants
  !-----------------------------------------------------------------------

  real(dp)     :: z_start_rec, z_end_rec, z_0, x_start_rec, x_end_rec, x_0

  integer(i4b)                           :: n_t, n1, n2        ! Number of x-values
  real(dp),    allocatable, dimension(:) :: x_t                ! Grid of relevant x-values
  real(dp),    allocatable, dimension(:) :: a_t                ! Grid of relevant a-values
  real(dp),    allocatable, dimension(:) :: eta_t              ! Grid of relevant eta-values

  integer(i4b)                           :: n_eta              ! Number of eta grid poins
  real(dp),    allocatable, dimension(:) :: x_eta              ! Grid points for eta
  real(dp),    allocatable, dimension(:) :: z_eta              ! Grid points for  z-values
  real(dp),    allocatable, dimension(:) :: a_eta              ! Grid points for a_eta-values
  real(dp),    allocatable, dimension(:) :: eta, eta2          ! Eta and eta'' at each grid point

  real(dp),    allocatable, dimension(:) :: rho_m              ! Matter density grid
  real(dp),    allocatable, dimension(:) :: rho_b              ! Baryon density grid
  real(dp),    allocatable, dimension(:) :: rho_r              ! Radiation density grid
  real(dp),    allocatable, dimension(:) :: rho_lambda         ! Dark energy density grid
  real(dp),    allocatable, dimension(:) :: rho_crit           ! Critical density grid
  real(dp),    allocatable, dimension(:) :: H_x                ! Hubble constant as function of x


contains

    !-----------------------------------------------------------------------
    ! Subroutines
    !-----------------------------------------------------------------------
  
  subroutine initialize_time_mod
    implicit none

    integer(i4b) :: i, n
    real(dp)     :: dx, dx_eta, x_eta1, x_eta2, a_init, a_end, eta_init, yp1, yp2, epsilon, steplength, steplength_min

    ! Define two epochs, 1) during and 2) after recombination.
  
    !-----------------------------------------------------------------------
    ! Initializing subroutine specific variables and constants
    !-----------------------------------------------------------------------

    n1          = 200                       ! Number of grid points during recombination
    n2          = 300                       ! Number of grid points after recombination
    n_t         = n1 + n2                   ! Total number of grid points

    z_start_rec = 1630.4d0                  ! Redshift of start of recombination
    z_end_rec   = 614.2d0                   ! Redshift of end of recombination
    z_0         = 0.d0                      ! Redshift today

    x_start_rec = -log(1.d0 + z_start_rec)  ! x of start of recombination
    x_end_rec   = -log(1.d0 + z_end_rec)    ! x of end of recombination
    x_0         = 0.d0                      ! x today
    
    n_eta       = 1000                      ! Number of eta grid points (for spline)
    a_init      = 1.d-10                    ! Start value of a for eta evaluation
    x_eta1      = log(a_init)               ! Start value of x for eta evaluation
    x_eta2      = 0.d0                      ! End value of x for eta evaluation

    yp1         = 1.d30                     ! natural spline "no derivatives at boundaries"
    yp2         = 1.d30
    epsilon     = 1.d-10
    ! Task: Fill in x and a grids


    !-----------------------------------------------------------------------
    ! Grids
    !-----------------------------------------------------------------------

    ! Allocating arrays
    allocate(x_t(n_t))
    allocate(a_t(n_t))
    allocate(x_eta(n_eta))
    allocate(a_eta(n_eta))
    allocate(z_eta(n_eta))
    allocate(eta(n_eta))
    allocate(eta2(n_eta))

    ! Filling in x-grid prior to recombination
    do i = 0, n1-1 
       x_t(i+1) = x_start_rec + i*(x_end_rec - x_start_rec)/(n1 - 1)
    end do

    ! Filling in x-grid post recombination
    do i = 1, n2
       x_t(n1+i) = x_end_rec + i*(x_0 - x_end_rec)/(n2)
    end do

    ! Creating a-grid using x = ln a
    a_t = exp(x_t)

    ! Task: 1) Compute the conformal time at each eta time step
    !       2) Spline the resulting function, using the provided "spline" routine in spline_1D_mod.f90
    
    dx_eta = (x_eta2 - x_eta1)/(n_eta - 1)              ! Step length of x_eta
    x_eta(1) = x_eta1                                   ! Initial condition of x_eta
    
    ! Filling in x_eta-grid
    do i = 1, n_eta-1
       x_eta(i+1) = x_eta(i) + dx_eta
    end do

    ! Computing a_eta and z_eta from x_eta
    a_eta = exp(x_eta)
    z_eta = 1.d0/a_eta - 1.d0


    !-----------------------------------------------------------------------
    ! Density Parameters
    !-----------------------------------------------------------------------

    ! Allocating the different arrays
    allocate(rho_m(n_eta))
    allocate(rho_b(n_eta))
    allocate(rho_r(n_eta))
    allocate(rho_lambda(n_eta))
    allocate(rho_crit(n_eta))

    ! Calculating the critical density
    do i = 1, n_eta
       rho_crit(i) = (3*get_H(x_eta(i))**2)/(8*pi*G_grav)
    end do

    rho_m = (Omega_m * rho_c) * a_eta**-3                ! Matter density
    rho_b = (Omega_b * rho_c) * a_eta**-3                ! Baryon density
    rho_r = (Omega_r * rho_c) * a_eta**-4                ! Radiation density
    rho_lambda = (Omega_lambda * rho_c)                  ! Dark energy density



    !-----------------------------------------------------------------------
    ! Eta Integration
    !-----------------------------------------------------------------------

    ! Computing the conformal time eta with the use of the ode_solver.f90 program
    eta(1)  = c*a_init/(H_0*sqrt(Omega_r))              ! Initial condition for eta
    steplength = abs((x_eta(1) - x_eta(2)))             ! Eta steplength
    steplength_min =  0                                 ! Minimum steplength

    ! Using the ODE-solver
    do i = 2, (n_eta)
       eta(i) = eta(i-1)
       call odeint(eta(i:i), x_eta(i-1), x_eta(i), epsilon, steplength, steplength_min, d_eta_da, bsstep, output)
    end do

    ! Forward Euler method (scrapped)
    ! do i=1, (n_eta-1)                                    
    !    eta(i+1) = eta(i) + (c/(get_H_p(x_eta(i+1))))*dx         
    ! end do
    
    ! Calling the spline function from spline_1D_mod.f90 program to do an interpolation
    call spline(x_eta, eta, yp1, yp2, eta2)


    !-----------------------------------------------------------------------
    ! Writing to file
    !-----------------------------------------------------------------------

    ! Writing quantities to 4 seperate files
    open(1, file = 'grids.dat', action = 'write', status = 'replace')
    open(2, file = 'densities.dat', action = 'write', status = 'replace')
    open(3, file = 'hubble.dat', action = 'write', status = 'replace')
    open(4, file = 'splined.dat', action = 'write', status = 'replace') 

    do i = 1, n_eta
       write(1,"(4(E17.8))") x_eta(i), a_eta(i), z_eta(i), eta(i)
       write(2,"(4(E17.8))") rho_m(i)/rho_crit(i), rho_b(i)/rho_crit(i), rho_r(i)/rho_crit(i), rho_lambda(i)/rho_crit(i)
       write(3,"(2(E17.8))") get_H(x_eta(i))
    end do

    do i = 1, n_t
       write(4,"(2(E17.8))") x_t(i), get_eta(x_t(i))
    end do

    close(1)
    close(2)
    close(3)
    close(4)

  end subroutine initialize_time_mod


  ! Subroutine for the eta derivative d(eta)/da  (may be smart to have it separate for future use)
  subroutine d_eta_da(x, eta, derivative)
    implicit none

    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: eta
    real(dp), dimension(:), intent(out) :: derivative
    
    derivative = c/get_H_p(x)

  end subroutine d_eta_da


  ! Subroutine for ODE-solver (output)
  subroutine output(x, y)
    implicit none

    real(dp),               intent(in) :: x
    real(dp), dimension(:), intent(in) :: y
    
  end subroutine output


  ! Task: Write a function that computes H at given x
  function get_H(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H
    real(dp)             :: a
    a = exp(x)

    get_H = H_0*sqrt((Omega_b + Omega_m)*a**-3.d0 + (Omega_r + Omega_nu)*a**-4.d0 + Omega_lambda)

  end function get_H


  ! Task: Write a function that computes H' = a*H  at given x
  function get_H_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H_p
    real(dp)             :: a
    a = exp(x)

    get_H_p = a*get_H(x)

  end function get_H_p


  ! Task: Write a function that computes dH'/dx at given x
  function get_dH_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dH_p
    real(dp)             :: a
    a = exp(x)
    get_dH_p = -((H_0**2)/(get_H_p(x))) * ((1/2)*(Omega_b + Omega_m)*a**-1 + Omega_r*a**-2 - Omega_lambda*a**2)
  end function get_dH_p


  ! Task: Write a function that computes eta(x), using the previously precomputed splined function
  function get_eta(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_eta
    real(dp)             :: a_in

    get_eta = splint(x_eta, eta, eta2, x_in)

  end function get_eta


end module time_mod
