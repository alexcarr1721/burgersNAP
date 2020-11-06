module functions
    use iso_fortran_env
    implicit none


    contains

    pure elemental function trap(x1,x2,f1,f2)
        ! Trapezoidal integration of function defined at 2 points ===============
        ! Input:
        !       x1      |   First x location
        !       x2      |   Second x location
        !       f1      |   Function value at x1
        !       f2      |   Function value at x2
        ! Output:
        !       trap    |   Trapezoidal approx of integral
        ! =======================================================================
        real(real64), intent(in) :: x1, x2, f1, f2
        real(real64) :: trap
        trap = ( f2 + f1 )*( x2 - x1 )/2.0
    end function trap

    subroutine trapz(f, x, start, err_code)
        ! Trapezoidal integration of function defined at 2 points ===============
        ! Input:
        !       f       |   Function to integrate
        !       x       |   Domain of function
        !       start   |   Initial Value of integral
        ! Output:
        !       f       |   Trapezoidal approx of integral (in-place)
        ! =======================================================================
        real(real64), intent(inout)             :: f(:)               
        real(real64), intent(in)                :: x(:)             
        integer(int64), intent(in), optional    :: start  
        integer(int64), intent(out), optional   :: err_code
        real(real64)                            :: initial_value
        real(real64), allocatable               :: f_old(:)
        integer(int64) :: i, j, k
        if ( present(err_code) ) err_code = 0

        if ( size(f) .ne. size(x) ) goto 100

        if ( present(start) ) then
          initial_value = start ! Set initial value of integral
        else
          initial_value = 0.0 ! Set to zero
        end if

        allocate( f_old(size(f)) )
        f_old = f

        do i = 1,size(f)
          if (i .eq. 1) then
            f(i) = initial_value
          else
            f(i) = trap(x(i-1),x(i),f_old(i-1),f_old(i)) + f(i-1)
          end if
        end do

        deallocate( f_old )
        goto 200

        100 continue
        if ( present(err_code) ) then
            err_code = 1
        else
            continue
        end if

        200 continue


    end subroutine trapz

    subroutine poisson(p, phi, theta, psi, dpsi, dsigma)
        ! =====================================================================
        ! Purpose: Compute poisson solution for Burgers equation              =
        !                                                                     =
        !                                                                     =
        !                                                                     =
        !         _________________________________________________________   =
        ! Inputs:|___Variable___|_______________Purpose____________________|  =
        !        |      p       |   Non-dimensional pressure               |  =
        !        |    theta     |   Non-dimensional time                   |  =
        !        |     psi      |                                          |  =
        !        |     dpsi     |                                          |  =
        !        |    dsigma    |                                          |  =
        !        |______________|__________________________________________|  =
        !          ________________________________________________________   =
        ! Outputs:|___Variable__|_______________Purpose____________________|  =
        !         |     p       |   Non-dimensional pressure               |  =
        !         |   theta     |   Atmospheric pressure at the ground     |  =
        !         |_____________|__________________________________________|  =
        !                                                                     =
        ! =====================================================================
        real(real64) :: p(:), phi(:), theta(:), psi(:), dpsi(:), dsigma
        ! Local variables =====================================================
        real(real64) :: theta_new(size(theta,dim=1))
        real(real64) :: phi_new(size(phi,dim=1))
        real(real64) :: dphi(size(phi,dim=1))
        integer(int64) :: i, j, k
        ! =====================================================================

        dphi = central_diff(phi, theta, "zero")

        ! Poisson solution ====================================================
        do i = 1,size(theta,dim=1)
          theta_new(i)    = theta(i) - dsigma*dphi(i)
          ! phi_new(i)      = phi(i) - 0.5*dsigma*( dphi(i)**2 )
        end do
        ! phi = phi_new
        theta = theta_new

        ! Extend solution and resample if beyond grid =========================
        ! =====================================================================


    end subroutine poisson

    subroutine branching(phi, theta)
        ! =====================================================================
        ! Purpose: Compute poisson solution for Burgers equation              =
        !                                                                     =
        ! Note:                                                               =
        !         _________________________________________________________   =
        ! Inputs:|___Variable___|_______________Purpose____________________|  =
        !        |      p       |   Non-dimensional pressure               |  =
        !        |    theta     |   Non-dimensional time                   |  =
        !        |     psi      |                                          |  =
        !        |     dpsi     |                                          |  =
        !        |    dsigma    |                                          |  =
        !        |______________|__________________________________________|  =
        !          ________________________________________________________   =
        ! Outputs:|___Variable__|_______________Purpose____________________|  =
        !         |     p       |   Non-dimensional pressure               |  =
        !         |   theta     |   Atmospheric pressure at the ground     |  =
        !         |_____________|__________________________________________|  =
        !                                                                     =
        ! =====================================================================
        real(real64)    :: phi(:), theta(:)
        ! Local variables =====================================================
        integer(int64)  :: branches, branch_size(size(phi,dim=1)+20)
        integer(int64)  :: branch_startend(2,size(phi,dim=1)+20) 
        real(real64)    :: phi_branches(size(phi,dim=1)+20,size(phi,dim=1)+20)
        real(real64)    :: theta_interp(size(theta,dim=1)+20)
        real(real64)    :: temp(size(phi,dim=1)+20)
        real(real64)    :: dtheta, start
        integer(int64)  :: indexMax
        integer(int64)  :: i, j, k
        real(real64)    :: theta_p(size(theta,dim=1)+20)
        real(real64)    :: phi_p(size(phi,dim=1)+20)
        ! =====================================================================

        dtheta = (maxval(theta,dim=1) - minval(theta,dim=1) )/( &
          1.0*(size(theta,dim=1) - 1) )
        do i = 1,size(theta_p,dim=1)
          if ( i .le. 10 ) then
            theta_p(i) = minval(theta,dim=1) - abs(i-11)*dtheta
            phi_p(i)   = 0.0
          else if ( i .gt. size(theta,dim=1)+10) then
            theta_p(i) = maxval(theta,dim=1) + (i - &
              (size(theta,dim=1)+10) )*dtheta
            phi_p(i)   = 0.0
          else 
            theta_p(i) = theta(i-10)
            phi_p(i)   = phi(i-10)
          end if
        end do

        ! Find branches in the solution ( if it is multivalued ) ==============
        branches = 0
        branch_size = 0
        branch_startend = 0
        branch_startend = 0
        do i = 1,size(theta_p,dim=1)-1
          if ( i .eq. 1) then
            ! Add the first branch ============================================
            branches = branches + 1
            branch_size(branches) = branch_size(branches) + 1
            branch_startend(1,branches) = i
            ! =================================================================
          else
            ! Add any subsequent branches if the solution is multivalued ======
            if ( theta_p(i+1) .ge. theta_p(i) ) then
              ! Theta is an increasing function in this interval ==============
              if ( i .eq. size(theta_p,dim=1) - 1 ) then
                ! End of theta vector =========================================
                if ( theta_p(i) .lt. theta_p(i-1) ) then
                  ! Was already on another branch, now there is a new branch at
                  ! the end. It starts at end-1 and goes to end. ==============
                  branches = branches + 1
                  branch_startend(1,branches) = i 
                  branch_size(branches) = branch_size(branches) + 2
                  branch_startend(2,branches) = i+1
                  ! ===========================================================
                else
                  ! Still on an increasing branch =============================
                  branch_size(branches) = branch_size(branches) + 2
                  branch_startend(2,branches) = i+1
                  ! ===========================================================
                end if
                ! =============================================================
              else if ( theta_p(i) .ge. theta_p(i-1) ) then
                ! Still on an increasing branch ===============================
                branch_size(branches) = branch_size(branches) + 1
                ! =============================================================
              else
                ! Transitioning from a decreasing branch ======================
                branches = branches + 1
                branch_size(branches) = branch_size(branches) + 1
                branch_startend(1,branches) = i
                ! =============================================================
              end if
              ! ===============================================================
            else
              ! Decreasing ====================================================
              if ( theta_p(i) .ge. theta_p(i-1) ) then
                ! Transitioning from an increasing branch =====================
                branch_size(branches) = branch_size(branches) + 1
                branch_startend(2,branches) = i
                ! =============================================================
              end if
              ! ===============================================================
            end if
          end if
        end do
        ! =====================================================================

        ! Assign values to individual branches ================================
        phi_branches = -9999999999999.9_real64 ! Initialize to really small number

        dtheta = (maxval(theta_p,dim=1) - minval(theta_p,dim=1) )/( &
          1.0*(size(theta_p,dim=1) - 1) )

        ! Create interpolation grid ===========================================
        start  = minval( theta_p, dim=1 )
        theta_interp = (/( (start + (i-1)*dtheta),i=1,size(theta_p,dim=1) )/)
        ! =====================================================================

        do i = 1,branches
          phi_branches(branch_startend(1,i):branch_startend(2,i),i) = &
            phi_p(branch_startend(1,i):branch_startend(2,i))
        end do

        ! Interpolate branches ================================================
        do i = 1,branches
          if ( branch_size(i) .eq. 1 ) then
            ! Re-interpolate onto finer grid? =================================
            write (*,*) "singleton"
          else
            ! Call interpolation subroutine ===================================
            call lininterp_dp(theta_p(branch_startend(1,i):&
              branch_startend(2,i)), &
              phi_branches(branch_startend(1,i):branch_startend(2,i),i), &
              theta_interp, temp, i, branches, -9999999999999.9_real64)
              phi_branches(:,i) = temp(:)
          end if
        end do
        ! =====================================================================

        ! Choose maximum ======================================================
        do i = 1,size(phi_p,dim=1)
          indexMax = maxloc( phi_branches(i,1:branches), dim=1 )
          phi_p(i)   = phi_branches(i,indexMax)
        end do
        ! =====================================================================

        dtheta = (maxval(theta,dim=1) - minval(theta,dim=1) )/( &
          1.0*(size(theta,dim=1) - 1) )
        start  = minval( theta, dim=1 )
        theta = (/( (start + (i-1)*dtheta),i=1,size(theta,dim=1) )/)
        call lininterp_dp(theta_interp, &
          phi_p, &
          theta, phi, 1_int64, 1_int64, -9999999999999.9_real64)

        ! ! Set theta and phi ===================================================
        ! theta = theta_interp(11:size(theta,dim=1)+10)
        ! phi   = phi_p(11:size(theta,dim=1)+10)
        ! ! =====================================================================

    end subroutine branching

    subroutine interp1d_branches(xData, yData, xVal, yVal, branch, branches,&
        out_of_bounds)
        ! Interpolation for Burger-Hayes ========================================
        ! Input:
        !       xData       |  a vector of the x-values of the data to be interpolated
        !       yData       |  a vector of the y-values of the data to be interpolated
        !       xVal        |  a vector of the x-values where interpolation should be performed
        ! Output:
        !       yVal          |  a vector of the resulting interpolated values
        ! =======================================================================
        real(real64), intent(in)    :: xData(:)
        real(real64), intent(in)    :: yData(:)
        real(real64), intent(in)    :: xVal(:)
        real(real64), intent(out)   :: yVal(:)
        integer(int64), intent(in)  :: branch, branches
        real(real64), intent(in)    :: out_of_bounds
        ! Local variables =====================================================
        integer(int64)  :: inputIndex, dataIndex
        integer(int64)  :: err_code
        real(real64)    :: minXdata, xRange, weight, maxXdata
        integer(int64)  :: i, j, k
        ! =====================================================================

        ! Possible checks on inputs could go here
        ! Things you may want to check:
        !   monotonically increasing xData
        !   size(xData) == size(yData)
        !   size(xVal) == size(yVal)
        if ( size(xData) .ne. size(yData) ) then
            write (*,*) "Error in interp1d():"
            write (*,*) "xData and yData must be the same length."
        end if
        if ( size(xVal) .ne. size(yVal) ) then
            write (*,*) "Error in interp1d():"
            write (*,*) "xVal and yVal must be the same length."
        end if
        do i = 1,size(xData)-1
            if ( ( xData(i+1) - xData(i) ) .lt. 0.0 ) then
                write (*,*) size(xData)
                write (*,*) "Error in interp1d():"
                write (*,*) "xData must be monotonically increasing."
            end if
        end do


        minXData = minval(xData) !xData(1)
        maxXData = maxval(xData) !xData(size(xData))
        xRange = maxXData - minXData

        ! write (*,*) minXdata, maxXdata

      do inputIndex = 1, size(xVal)
        ! Set out of range values to 0
        if (xVal(inputIndex) .lt. minXData) then
          if ( branch .eq. 1 ) then
            yVal(inputIndex) = 0.0
          else
            yVal(inputIndex) = out_of_bounds
          end if
        else if (xVal(inputIndex) .gt. maxXData) then
          if ( branch .eq. branches ) then
            yVal(inputIndex) = 0.0
          else
            yVal(inputIndex) = out_of_bounds
          end if
        else

          ! If dataIndex .eq. 1, check xData(inputIndex +- 1)
          if ( inputIndex .eq. size(xVal) ) goto 200
          if ( inputIndex .eq. 1 ) goto 200
          if ( inputIndex .ge. size(xData)-1 ) goto 200
          if ( (xData(inputIndex+1) .ge. xVal(inputIndex)) .and. &
            (xData(inputIndex) .le. xVal(inputIndex)) ) then
            ! Skip to interpolation
            dataIndex = inputIndex
            goto 100
          elseif ( (xData(inputIndex) .ge. xVal(inputIndex)) .and. &
            (xData(inputIndex-1) .le. xVal(inputIndex)) ) then
            ! Skip to interpolation
            dataIndex = inputIndex
            goto 100
          else
            ! Continue
            continue
          end if

          200 continue
          ! this will work if x is uniformly spaced, otherwise increment
          ! dataIndex until xData(dataIndex+1)>xVal(inputIndex)
          dataIndex = 0
          findloop: do i = 1, size(xData)-1
            dataIndex = dataIndex + 1
            ! If dataIndex .eq. 1, check xData(inputIndex +- 1)
            if ( (xData(dataIndex+1) .ge. xVal(inputIndex)) .and. &
              (xData(dataIndex) .le. xVal(inputIndex)) ) then
              ! Exit loop
              go to 100
            else
              ! Continue loop
              continue
            end if
          end do findloop

          100 continue
          if ( yVal(dataIndex+1) .lt. -5e5 ) write (*,*) "Yikes"
          ! write (*,*) yData(dataIndex)
          if ( xData(dataIndex+1) .eq. xData(dataIndex) ) then ! xVal must be equal to xData
            yVal(inputIndex) = yData(dataIndex)
          else
            weight = (xVal(inputIndex) - &
              xData(dataIndex))/(xData(dataIndex+1)-xData(dataIndex))
            yVal(inputIndex) = (1.0-weight)*yData(dataIndex) + &
              weight*yData(dataIndex+1)
          end if
        end if
      end do

    end subroutine interp1d_branches

    function central_diff(f, x, bc)
        ! Central difference of a function f with domain x ======================
        ! Input:
        !       f       |   function array
        !       x       |   domain
        !       bc      |   boundary condition (optional)
        ! Output:
        !       central_diff    |   Central diff approx of derivative df/dx
        ! =======================================================================
        real(real64), intent(in)                   :: f(:)
        real(real64), intent(in)                   :: x(:)
        character(len=*), intent(in), optional :: bc
        real(real64)                               :: central_diff(size(f))
        integer(int64)                              :: i, j, k
        if ( present(bc) ) then
            if ( trim(bc) .eq. 'periodic') then
                goto 100
            else if ( trim(bc) .eq. 'zero') then
                goto 200
            else if ( trim(bc) .eq. 'derivative' ) then
                goto 300
            else
                goto 300
            end if
        else
            goto 100
        end if

        100 continue ! Periodic boundary
        do i = 1, size(x)
          if ( i .eq. 1 ) then
            central_diff(i) = ( f(i+1) - f(size(x)) )/( 2.0*( x(i+1) - x(i) ) )
          else if ( i .eq. size(x) ) then
            central_diff(i) = ( f(i-1) - f(1) )/( 2.0*(x(i) - x(i-1)) )
          else
            central_diff(i) = ( f(i+1) - f(i-1) )/( 2.0*(x(i+1) - x(i) ) )
          end if
        end do
        goto 400

        200 continue
        do i = 1, size(x)
          if ( i .eq. 1 ) then
            central_diff(i) = ( f(i+1) - f(i+1) )/( 2.0*( x(i+1) - x(i) ) )
          else if ( i .eq. size(x) ) then
            central_diff(i) = ( f(i-1) - f(i-1) )/( 2.0*(x(i) - x(i-1)) )
          else
            central_diff(i) = ( f(i+1) - f(i-1) )/( 2.0*(x(i+1) - x(i) ) )
          end if
        end do
        goto 400

        300 continue
        do i = 1, size(x)
          if ( i .eq. 1 ) then
            central_diff(i) = ( f(i+1) - f(i) )/( ( x(i+1) - x(i) ) )
          else if ( i .eq. size(x) ) then
            central_diff(i) = ( f(i) - f(i-1) )/( (x(i) - x(i-1)) )
          else
            central_diff(i) = ( f(i+1) - f(i-1) )/( 2.0*(x(i+1) - x(i) ) )
          end if
        end do
        goto 400

        400 continue

    end function central_diff

    subroutine search_table(x,xData)
      real(real64) :: x, xData(:)

      

    end subroutine search_table

  !   subroutine hunt(x, xData, jlo)
  !     integer(idp), intent(inout) :: jlo
  !     real(dp), intent(in)        :: x
  !     real(dp), intent(in)        :: xData(:)
  !     ! Local variables =====================================================
  !     integer(idp)                :: n, inc, jhi, jm
  !     logical                     :: ascnd
  !     ! =====================================================================

  !     n = size(xData,dim=1)
  !     ascnd = ( xData(n) .ge. xData(1) )

  !     if ( (jlo .le. 0) .or. (jlo .gt. n) ) then
  !         jlo = 0
  !         jhi = n + 1
  !     else
  !         inc = 1
  !         if ( (x .ge. xData(jlo) ) .eqv. ascnd ) then
  !             do
  !                 jhi = jlo + inc
  !                 if (jhi .gt. n) then
  !                     jhi = n + 1
  !                     exit    
  !                 else
  !                     if ( (x .lt. xData(jhi) ) .eqv. ascnd ) exit 
  !                     jlo = jhi
  !                     inc = inc + inc
  !                 end if
  !             end do
  !         else
  !             jhi = jlo 
  !             do  
  !                 jlo = jhi - inc 
  !                 if ( jlo .lt. 1 ) then
  !                     jlo = 0
  !                     exit 
  !                 else
  !                     if ( (x .gt. xData(jlo) ) .eqv. ascnd ) exit
  !                     jhi = jlo 
  !                     inc = inc + inc 
  !                 end if 
  !             end do
  !         end if
  !     end if

  !     ! Bisection
  !     do  
  !         if ( jhi-jlo .le. 1 ) then 
  !             if ( x .eq. xData(n) ) jlo = n-1
  !             if ( x .eq. xData(1) ) jlo = 1
  !             exit 
  !         else
  !             jm = (jhi + jlo)/2
  !             if ( (x .ge. xData(jm)) .eqv. ascnd ) then
  !                 jlo = jm
  !             else
  !                 jhi = jm 
  !             end if
  !         end if 
  !     end do

  ! end subroutine hunt

  function locate(x,xData)
      ! =====================================================================
      !   Purpose:    Returns a value (locate) such that x is between       =
      !               xData(locate) and xData(locate + 1).                  =
      !                                                                     =
      ! =====================================================================
      real(real64), intent(in)    :: x
      real(real64), intent(in)    :: xData(:)
      integer(int64)            :: locate
      ! Local variables =====================================================
      integer(int64)            :: n, jl, jm, ju
      logical                 :: ascnd 
      ! =====================================================================

      n = size(xData,dim=1)
      ascnd = ( xData(n) .ge. xData(1) )
      jl = 0                      ! Lower limit
      ju = n+1                    ! Upper limit
      do while ( ju - jl > 1 )    ! Bisection
          jm = (ju + jl)/2        ! Compute the midpoint
          if ( ascnd .eqv. (x .ge. xData(jm) ) ) then
              jl = jm             ! Replace lower bound
          else
              ju = jm             ! Replace upper bound
          end if
      end do

      if ( x .eq. xData(1) ) then
          locate = 1
      else if ( x .eq. xData(n) ) then
          locate = n-1
      else
          locate = jl
      end if

  end function locate

  subroutine lininterp_dp(xData, yData, x, y, n, m, out_of_bounds)
      real(real64), intent(in) :: xData(:)
      real(real64), intent(in) :: yData(:)
      real(real64), intent(in) :: x(:)
      real(real64), intent(out):: y(:)
      integer(int64), intent(in) :: n ! Branch number
      integer(int64), intent(in) :: m ! number of branches
      real(real64), optional   :: out_of_bounds 
      ! Local
      real(real64)             :: oob, xMin, xMax
      integer(int64)         :: i, j, k, jl

      if ( present(out_of_bounds) ) then
          oob = out_of_bounds
      else    
          oob = -9999999999999.9_real64
      end if

      ! Insert checks for monotonicity of xData


      xMin = minval( xData )
      xMax = maxval( xData )
      ! Loop through x
      do i = 1,size(x,dim=1)
          if ( x(i) .lt. xMin ) then
              if ( n .eq. 1 ) then
                  y(i) = 0.0_real64
                  go to 100
              else
                  y(i) = oob 
                  go to 100
              end if
          else if ( x(i) .gt. xMax ) then
              if ( n .eq. m ) then 
                  y(i) = 0.0_real64
              else
                  y(i) = oob 
              end if
          else
              ! Search table to find location of x in xData domain
              jl = locate(x(i), xData)
              ! Interpolate?
              if ( yData(jl) .le. oob ) write (*,*) "Out of bounds y value."
              if ( yData(jl+1) .le. oob ) write (*,*) "Out of bounds y value."
              y(i) = yData(jl) + ( (x(i) - xData(jl) )/(xData(jl+1) &
                  - xData(jl)) )*(yData(jl+1) - yData(jl))
              
          end if
          100 continue
      end do

  end subroutine lininterp_dp

  subroutine waveform_parameters(p, t, c0, rho0, gamma, pa, omega0, xbar)
    !**********************************************************************
    ! Purpose:                                                            *
    !           Obtain parameters of the input waveform like the period,  *
    !           Amplitude, etc.                                           *         
    !                                                                     *  
    ! Inputs:    _________________________________________________________* 
    !           |___Variable___|_______________Purpose____________________|  
    !           |       p      |   Pressure waveform                      |  
    !           |       t      |   Time                                   | 
    !           |      c0      |   Speed of sound                         | 
    !           |    rho0      |   Density                                |
    !           |   gamma      |   Ratio of specific heats                |
    !           |______________|__________________________________________|  
    !                                                                     *
    ! Outputs:   _________________________________________________________*
    !           |___Variable____|_______________Purpose___________________|  
    !           |      pa       |   Maximum of the absolute value of p    |
    !           |  omega0       |   2pi/period of signal                  |
    !           |    xbar       |   Shock formation distance for equi-    |
    !           |               |   valent sinusoidal signal              |
    !           |_______________|_________________________________________|  
    !                                                                     *
    !**********************************************************************
    implicit none
    real(real64), intent(in)    :: p(:), t(:), c0, rho0, gamma
    real(real64), intent(out)   :: pa, omega0, xbar
    ! Local variable ******************************************************
    integer(int64)              :: i, last, start
    real(real64)                :: period, eps, p_reverse(size(p,dim=1))
    real(real64)                :: f(size(p,dim=1))
    real(real64)                :: beta
    ! Parameters **********************************************************
    real(real64), parameter     :: pi = 3.141592653589793
    !**********************************************************************                

    ! Maximum *************************************************************
    pa = maxval( abs(p) )
    !**********************************************************************

    ! Store reverse array *************************************************
    do i = 1,size(p,dim=1)
        p_reverse(i) = p(size(p,dim=1) - i + 1)
    end do
    !**********************************************************************

    ! Find period *********************************************************
    eps = 0.01                                                         ! 1%
    do i = 1,size(p,dim=1)
        if ( abs(p(i)) .ge. eps*pa ) then
            start = i
            go to 200
        else
            continue
        end if
    end do

    200 continue
    do i = 1,size(p,dim=1)
        if ( abs(p_reverse(i)) .ge. eps*pa ) then
            last = size(p,dim=1) - i + 1
            go to 300
        else
            continue
        end if
    end do
    !**********************************************************************

    300 continue

    ! omega0 **************************************************************
    f = central_diff(p, t, bc="periodic")
    period = t(last) - t(start)
    omega0 = 2.0*pi/period
    !**********************************************************************

    ! Consider the wave to have begun initially as a sine wave, compute ***
    ! the shock formation distance using this assumption ******************
    beta = 0.5*(gamma + 1.0)
    xbar = rho0*(c0**3)/(beta*maxval( abs(f) ))
    !**********************************************************************

  end subroutine waveform_parameters

end module functions