module cl_mod
  use healpix_types
  use evolution_mod
  use sphbess_mod
  implicit none

  ! ======================
  ! Global parameters
  ! ======================

  real(dp),     pointer,     dimension(:,:)     :: j_l, j_l2, S, S2
  real(dp),     pointer,     dimension(:)       :: x_highres, k_highres, l_highres, cl_highres
  real(dp),     allocatable, dimension(:)       :: z_spline, j_l_spline, j_l_spline2, callin_int


contains

  ! ======================
  ! Subroutines
  ! ======================

  subroutine compute_cls
    implicit none

    ! ------------------------------------------------------------
    ! Function which computes the CMB power spectrum cl
    ! ------------------------------------------------------------

    ! Parameters
    integer(i4b)       :: i, j, k, l, l_num, l_max, n_spline
    real(dp)           :: dx, dk,  S_func, j_func, z, x_min, x_max, k_min, k_max, n_spec
    real(dp)           :: z_start, z_stop, yp1, yp2, integral

    logical(lgt)       :: exist
    character(len=128) :: filename

    ! Arrays
    integer(i4b), allocatable, dimension(:)       :: ls
    real(dp),     allocatable, dimension(:,:)     :: Integrand_Theta, Integrand_cl
    real(dp),     pointer,     dimension(:)       :: cls, cls2, ls_dp, cl_int
    real(dp),     pointer,     dimension(:,:,:,:) :: S_coeff
    real(dp),     allocatable, dimension(:,:)     :: Theta_l

    ! Spline variables
    yp1     = 1.d30
    yp2     = 1.d30
    
    ! Length of spline array
    n_spline = 5400

    l_num = 44
    l_max = 1200


    ! ---------------------------------------------------------
    ! Allocating arrays to be used during calculations (many)
    ! ---------------------------------------------------------

    ! Arrays which come from get_higres_source_function
    allocate(x_highres(n_x_highres))
    allocate(k_highres(n_k_highres))
    allocate(S(n_x_highres, n_k_highres))
    
    ! Array containing the different l values to be used
    allocate(ls(l_num))

    ! High resolution and double precsision arrays
    allocate(ls_dp(l_num))
    allocate(l_highres(l_max))
    allocate(cl_highres(l_max))

    ! Spline related arrays (Bessel)
    allocate(z_spline(n_spline))    ! Note: z is *not* redshift, but simply the dummy argument of j_l(z)
    allocate(j_l(n_spline, l_num))
    allocate(j_l2(n_spline, l_num))

    ! Callin figure comparison check
    allocate(callin_int(n_x_highres))

    ! Integration related, Transfer function and its integrand 
    allocate(Theta_l(l_num, n_k_highres))
    allocate(Integrand_Theta(l_num, n_x_highres))
    
    ! Integration related, cl function and its integrand
    allocate(cl_int(l_num))
    allocate(cls(l_num))
    allocate(cls2(l_num))
    allocate(Integrand_cl(l_num, n_k_highres))

    ! Opening files
    !open(29, file = 'bessel.unf', form = 'unformatted', action = 'write', status = 'replace')
    open(30, file = 'cl_best.dat', action = 'write', status = 'replace')
    !open(31, file = 'source.unf', form = 'unformatted', action = 'write', status = 'replace')
    !open(32, file = 'transfer.unf', form = 'unformatted', action = 'write', status = 'replace')
    !open(33, file = 'callin.dat', action = 'write', status = 'replace')
    !-------------------------
    ! Begining calculations
    !-------------------------

    ! Set up which l's to compute
    ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
         & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
         & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)


    ! Converting ls from int to double (to be used at the end) and making highres arrays
    do l = 1, l_num

       ls_dp(l) = ls(l)

    end do

    ! Making double precission array from 1 to l_max
    do l = 1, l_max

       l_highres(l) = l

    end do

    ! Extracting the highresolution source function from evo_mod 
    call get_hires_source_function(x_highres, k_highres, S)

    ! Initializing z array to be splined going from "0" to 3500 (not including 0)
    z_start  = 0.d0
    z_stop   = 3500.d0 
    do i = 1, n_spline

       z_spline(i) = z_start + (i-1.d0)*((z_stop - z_start)/(n_spline-1.d0))
       
    end do
    
    ! Precoming the bessel functions
    do i = 1, n_spline

       do l = 1, l_num
          if (z_spline(i) > 0.d0) then
             call sphbes(ls(l), z_spline(i), j_l(i,l))
          end if
       end do

    end do

    ! Splining the bessel functions
    do l = 1, l_num

       call spline(z_spline, j_l(:,l), yp1, yp2, j_l2(:,l))

    end do
    
    ! Writing bessel functions to file
    
!    write(29) j_l
!    write(29) j_l2

    ! Computing integrand plot similar to callins
!    j = locate_dp(k_highres, 340.d0*H_0/c)
 
!    do i = 1, n_x_highres

 !      callin_int(i) = S(i,j)*j_lfunction(17, x_highres(i), k_highres(j))
 !      write(33, '(1(E17.8E3))') callin_int(i)

!    end do



    ! -------------------------------------------------
    ! Computing the Transfer function and the cl's
    ! -------------------------------------------------
    
    ! Computing variables for integration below (dx and dk)
    x_min = x_highres(1)
    x_max = x_highres(n_x_highres)
    k_min = k_highres(1)
    k_max = k_highres(n_k_highres)

    dx = (x_max - x_min)/n_x_highres
    dk = (k_max - k_min)/n_k_highres

    ! Loop of l's
    do l = 1, l_num

       write(*,*) l

       ! Initializing the transfer function and cls as 0 to avoid random memory usage
       Theta_l(l,:) = 0.d0
       cl_int(l)    = 0.d0

       ! Loop over k's
       do k = 1, n_k_highres
          
          ! Loop over x's
          do i = 1, n_x_highres

             ! Computing the integrand in the transfer function Theta_l
             Integrand_Theta(l,i) = S(i,k)*j_lfunction(l,x_highres(i), k_highres(k))
             
          end do
          
          do i = 1, n_x_highres
             
             ! Summing the integrands resulting in unintegrated Theta
             Theta_l(l,k) = Theta_l(l,k) + Integrand_Theta(l,i)

           end do

          ! "Integrating" up Theta_l
          Theta_l(l,k) = Theta_l(l,k)*dx

          !Computing P(k) * (Theta_l^2/k) by first computing the integrand of cl 
          Integrand_cl(l,k) = (((c*k_highres(k))/H_0)**(n_s - 1.d0))*((Theta_l(l,k)**(2.d0))/k_highres(k))
          cl_int(l) = cl_int(l) + Integrand_cl(l,k)

       end do
       
       ! Task: Integrate P(k) * (Theta_l^2 / k) over k to find un-normalized C_l's

       ! Task: Store C_l in an array. Optionally output to file
       cls(l) = cl_int(l)*((ls(l)*(ls(l) + 1.d0))/(2.d0*pi))*dk

    end do
    

    ! Task: Spline C_l's found above, and output smooth C_l curve for each integer l

    ! Splining the resultin cl's found above to give us a smooth c_l curve for wanted l's
    call spline(ls_dp, cls, yp1, yp2, cls2)

    ! Writing to file
    do l = 1, ls(l_num)

       cl_highres(l) = splint(ls_dp, cls, cls2, l_highres(l))

       write(30, '(2(E17.8))') l_highres(l), cl_highres(l)

    end do

 !   write(31) x_highres
 !   write(31) k_highres
 !   write(31) S
 !   write(32) Integrand_Theta
 !   write(32) Theta_l 
 !   write(32) Integrand_cl
 !   write(32) cl_int

    ! Closing files
 !   close(29)
    close(30)
 !   close(31)
 !   close(32)
 !   close(33)

  end subroutine compute_cls
  
  function j_lfunction(l, x, k)
    implicit none
    integer(i4b),     intent(in) :: l
    real(dp),         intent(in) :: x, k
    real(dp)                     :: j_lfunction
    
    j_lfunction = splint(z_spline, j_l(:,l), j_l2(:,l), k*(get_eta(0.d0) - get_eta(x)))

  end function j_lfunction

end module cl_mod
