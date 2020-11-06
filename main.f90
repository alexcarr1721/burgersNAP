! =============================================================================
!               NAP: Nonlinear Acoustic Propagation code                      |
!                                                                             |
!   Author: Alexander N. Carr                                                 |
!   Institution: University of Florida                                        |
!   Occupation: Ph.D. student in the Mechanical and Aerospace Engineering Dpt.|
!   Group: Theoretical Fluid Dynamics and Turbulence Group                    |
!   Contact: alexcarr.1721@gmail.com (Please include "NAP code" in subject)   |
!                                                                             |
!   Support: This work is supported by NASA Fellowship Grant No. <insert>     |
! =============================================================================

program main 
    use iso_fortran_env
    use functions
    implicit none 
    ! Variables ===============================================================
    real(real64), allocatable :: p(:), x(:), tau(:), pnond(:)
    real(real64), allocatable :: phi(:), psi(:), dpsi(:), theta(:), sigma(:)
    real(real64)  :: c0, rho0, p0, gamma
    ! real(real64) :: p(1024), x(1024), tau(1024), pnond(1024)
    ! real(real64), parameter :: c0     = 340.0
    ! real(real64), parameter :: rho0   = 1.25
    ! real(real64), parameter :: p0     = 101325.0
    ! real(real64), parameter :: gamma  = 1.4
    ! real(real64), parameter :: f      = 1.0
    real(real64), parameter :: pi     = 3.141592653589793
    ! real(real64), parameter :: omega0 = f*pi*2
    ! Non-dimensional parameters ==============================================
    real(real64)            :: omega0, beta, pa, x_bar
    real(real64)            :: eps
    real(real64)            :: k0, start, fmax, dx
    integer(int64)          :: i, j, k
    ! File read variables =====================================================
    character(len=80)       :: trash1
    integer(int64)          :: xsize, tsize
    real(real64)            :: propagation_distance
    ! =========================================================================

    ! Read file ===============================================================
    open (unit = 10, file = "domain.dat")
    read (10,*) trash1, tsize
    read (10,*) trash1, xsize
    read (10,*) trash1, propagation_distance
    read (10,*) trash1, c0
    read (10,*) trash1, rho0
    read (10,*) trash1, p0
    read (10,*) trash1, gamma
    close(10)
    ! =========================================================================

    ! Allocate ================================================================
    allocate ( p(tsize), x(xsize), tau(tsize), pnond(tsize) )
    allocate ( phi(tsize), psi(tsize), dpsi(tsize), theta(tsize), sigma(xsize))
    ! =========================================================================

    ! Read waveform ===========================================================
    open (unit = 10, file = "waveform.dat")
    do i = 1,tsize
      read (10,*) tau(i), p(i)
    end do
    close(10)
    ! =========================================================================

    ! Non-dimensional parameters ==============================================
    call waveform_parameters(p, tau, c0, rho0, gamma, pa, omega0, x_bar)
    ! =========================================================================

    ! X =======================================================================
    dx = propagation_distance*x_bar/xsize
    x = (/ ( (i-1)*( dx/(xsize - 1) ), i=1,xsize  ) /)


    ! Domain ==================================================================
    ! eps = p0/(rho0*(c0**2))
    ! beta = 0.5*(gamma+1)
    ! k0 = (2.0*pi*f)/c0
    ! x_bar = 1.0/(beta*eps*k0)
    ! write (*,*) x_bar
    ! start  = -pi/(k0*c0)
    ! tau = (/ ( ( start + (i-1)*( (pi/(k0*c0) -start)/(size(tau,dim=1) -1) ) &
    !   ),i=1,size(tau,dim=1) ) /)

    ! ! Initial Condition =======================================================
    ! pa = p0
    ! ! do i = 1,size(tau,dim=1)
    ! !     p(i) = pa*sin(omega0*tau(i))
    ! ! end do
    ! do i = 1,size(tau,dim=1)
    !   ! Construct 1st section of N-wave
    !   if ( i .le. size(tau,dim=1)/2 ) then
    !     p(i) = -1.0*(i-1)*1.0/(1.0*(size(tau,dim=1)/2) )*pa
    !   else
    !     p(i) = pa*(size(tau,dim=1) - i)*1.0/(1.0*(size(tau,dim=1)/2 + 1) )
    !   end if
    ! end do
    ! pnond = p/pa
    ! ! =========================================================================
    ! ! fmax = maxval( abs(central_diff(p, tau, "derivative")) )
    ! ! x_bar = rho0*(c0**3)/(beta*fmax)
    ! write(*,*) x_bar
    ! start  = 0.0
    ! x = (/ ( ( start + (i-1)*( (10.0*x_bar - start)/(size(x,dim=1) - 1) ) &
    !   ),i=1,size(x,dim=1) ) /)
    ! ! start  = -pi/(k0*c0)
    ! ! tau = (/ ( ( start + (i-1)*( (pi/(k0*c0) -start)/(size(tau,dim=1) -1) ) &
    ! !   ),i=1,size(tau,dim=1) ) /)
    ! sigma = x/x_bar 
    ! theta = omega0*tau

    ! write (*,*) ( sigma(2)-sigma(1) )/( theta(2) - theta(1) )
    ! ! =========================================================================

    ! ! Initial Condition =======================================================
    ! pa = p0
    ! ! do i = 1,size(tau,dim=1)
    ! !     p(i) = pa*sin(omega0*tau(i))
    ! ! end do
    ! do i = 1,size(tau,dim=1)
    !   ! Construct 1st section of N-wave
    !   if ( i .le. size(tau,dim=1)/2 ) then
    !     p(i) = -1.0*(i-1)*1.0/(1.0*(size(tau,dim=1)/2) )*pa
    !   else
    !     p(i) = pa*(size(tau,dim=1) - i)*1.0/(1.0*(size(tau,dim=1)/2 + 1) )
    !   end if
    ! end do
    ! pnond = p/pa
    ! ! =========================================================================

    ! Non-dimensional forms ===================================================
    pnond = p!p/pa 
    sigma = x!x/x_bar 
    theta = tau!omega0*tau
    ! =========================================================================

    ! Compute psi and dpsi ====================================================
    do i = 1,size(tau,dim=1)
        psi(i) = 0.5*(pnond(i)**2)
        dpsi(i) = pnond(i)
    end do
    ! =========================================================================

    ! Integrate to get potential ==============================================
    call trapz(pnond, theta)
    phi = pnond
    ! call trapz(p, tau)
    ! phi = p*omega0/pa
    ! =========================================================================

    ! Progress forward in space ===============================================
    do i = 2,size(sigma,dim=1)

      ! write (*,*) minval( phi )
      ! Compute poisson solution ==============================================
      call poisson(pnond, phi, theta, psi, dpsi, sigma(i)-sigma(i-1))
      ! =======================================================================
      ! write (*,*) minval( phi )

      ! Find increasing branches and interpolate ==============================
      call branching(phi, theta)
      ! =======================================================================
      ! write (*,*) minval( phi )

      ! Compute Pressure ======================================================
      pnond = central_diff(phi, theta, "zero")
      ! =======================================================================

      ! Compute psi and dpsi ==================================================
      do j = 1,size(theta,dim=1)
        psi(j) = 0.5*(pnond(j)**2)
        dpsi(j) = pnond(j)
      end do
      ! =======================================================================

      write (*,*) maxval( pnond ), pnond(1), pnond(1024), x(i)

    end do
    ! =========================================================================

    ! Write to file ===========================================================
    open(10, file = 'out.dat', status = 'old')  
    do i=1,size(theta,dim=1)  
      write(10,*) theta(i), pnond(i)   
    end do

end program main