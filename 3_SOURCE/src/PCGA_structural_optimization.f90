    module PCGA_struct_param_optimization

    !*********************************************************************************************************************************************************
    ! Method for optimizing the structural parameters. The PCGA parts should have been able to implement in the original structural optimization, however,
    ! it was not possible to make it work. Therefore, the original structural optimization was copied and changed such it fit with the PCGA optimizations.
    !*********************************************************************************************************************************************************
    
    use bayes_pest_control
    use model_input_output
    use bayes_pest_finalize
    use error_message
    use utilities
    use make_kernels
    use bayes_matrix_operations
    use jupiter_input_data_support
    use nelder_mead_simplex_routines
    use PCGA_solve_linear_equations
    use PCGA_low_rank_app
    use PCGA_param_trans

    contains

    subroutine PCGA_marginal_struct_param_optim (d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,errstruc,miostruc,cv_MIO, d_MIO,d_MOD,d_ANI,it_bga,n)

    implicit none
    type(kernel_XQR),     intent(inout)     :: d_XQR
    type(cv_observ),      intent(inout)     :: cv_OBS
    type (cv_minout),     intent(in)        :: cv_MIO
    type(d_observ),       intent(inout)     :: d_OBS
    type(d_algorithmic),  intent(inout)     :: d_A
    type(cv_algorithmic), intent(inout)     :: cv_A
    type(d_prior_mean),   intent(in)        :: d_PM
    type(cv_param),       intent(inout)     :: cv_PAR
    type(d_param),        intent(inout)     :: d_PAR
    type(cv_struct),      intent(inout)     :: cv_S
    type(d_struct),       intent(inout)     :: d_S
    type(cv_prior_mean),  intent(in)        :: cv_PM
    character (len=ERRORWIDTH)              :: retmsg

    !Introduced for the PCGA part
    type (err_failure_struc) , intent(inout)            :: errstruc
    type (mio_struc)         , intent(inout)            :: miostruc
    type (d_comlin)          , intent(inout)            :: d_MOD
    type(d_anisotropy)       , intent(in)               :: d_ANI
    type (d_minout)          , intent(in)               :: d_MIO
    type(Q0_compr)           , intent(inout), pointer   :: Q0_All(:)

    integer              :: it_bga, i, j, k, p, imax
    integer              :: nQ0 = 0  !Dimension of Q0_All(:) [0] in case of no compression [cv_PAR%p] in case of compression
    double precision     :: stc = 1. !Parameter that control values of step. step = stc*initial values
    integer ( kind = 4 ) :: n        !Indicates the number of pars to estimate.
    integer ( kind = 4 ) :: icount   !Number of function evaluation used (output)
    integer ( kind = 4 ) :: ifault   !Error indicator (output)
    integer ( kind = 4 ) :: konvge
    integer ( kind = 4 ) :: numres   !Number of restarts (output)
    real    ( kind = 8 ) :: reqmin
    real    ( kind = 8 ) :: start(n)
    real    ( kind = 8 ) :: step(n)
    real    ( kind = 8 ) :: xmin(n)
    real    ( kind = 8 ) :: ynewlo

    !***********************************************************************************************************************************
    !Start is the starting point for the iteration, and must contain the power transformed values if power transformation is required.
    !Step determine size and shape of initial simplex and must reflect the units of the variables. step = stc * initial values.
    !In case of power transformation step is calculated to be the same of the not tansformed case.
    !***********************************************************************************************************************************

    start = d_S%struct_par_opt_vec
    step  = stc * d_S%struct_par_opt_vec

    !Loop to power transform where necessary
    if ((maxval(cv_S%trans_theta).eq.1).or.(d_S%trans_sig.eq.1)) then !This means we have at least a theta or sig to power transform
        k = 0
        if ((maxval(cv_S%struct_par_opt).eq.1)) then
            do i = 1, cv_PAR%p
                if (cv_S%struct_par_opt(i).eq.1)  then
                    do j = 1,cv_S%num_theta_type (i)
                        k = k + 1
                        if (cv_S%trans_theta(i).eq.1) then
                            start(k) = cv_S%alpha_trans(i)*((d_S%struct_par_opt_vec(k)**(1./cv_S%alpha_trans(i)))-1) !Forward-trans
                            step(k) = cv_S%alpha_trans(i)*((((stc+1)*d_S%struct_par_opt_vec(k))**(1./cv_S%alpha_trans(i)))-1) - &
                                cv_S%alpha_trans(i)*((d_S%struct_par_opt_vec(k)**(1./cv_S%alpha_trans(i)))-1)
                        endif
                    enddo
                endif
            enddo
        endif

        if ((d_S%sig_opt.eq.1).and.(d_S%trans_sig.eq.1)) then
            k = k + 1
            start(k) = d_S%alpha_trans_sig *((d_S%struct_par_opt_vec(k)**(1./d_S%alpha_trans_sig))-1) !Forward-trans
            step(k) = d_S%alpha_trans_sig *((((stc+1)*d_S%struct_par_opt_vec(k))**(1./d_S%alpha_trans_sig))-1) - &
                d_S%alpha_trans_sig *((d_S%struct_par_opt_vec(k)**(1./d_S%alpha_trans_sig))-1)
        endif
    endif

    konvge = 1
    reqmin = 1.0D-2 !Need to be tested
    if (cv_A%Q_compression_flag.ne.0) nQ0 = cv_PAR%p

    !Using downhill simplex / Nelder–Mead method to minimize ynewlo
    call PCGA_nelmin_sp(n,start,xmin,ynewlo,reqmin,step,konvge,cv_A%it_max_structural,icount,numres,ifault, &
        & d_XQR, Q0_all,cv_OBS, d_OBS, cv_A, d_A, d_PAR, cv_S, d_S, d_PM, cv_PAR,cv_PM,nQ0, errstruc, miostruc, cv_MIO, d_ANI, d_MIO, d_MOD)


    if (ifault.eq.2) then !Maximum number of iterations has exceeded --> Warning
        write(retmsg,10) it_bga
10      format('Warning: Maximum number of iterations exceeded in structural parameter optimization procedure during bgaPEST iteration',i4, &
            & '. The vector that gives the minimum obj. funct. during the procedure will be considered.')
        call utl_writmess(6,retmsg)
    endif

    cv_S%str_obj_fun = ynewlo ! Minimum value of the objective function

    !****************************************************************************************************************************************
    !------ Here we reform d_S%theta and d_S%sig overwriting the elements optimized by Nelder-Mead and leaving unchanged the others. --------
    !---------- The values that minimize the structural parameter objective function are also assigned to d_S%struct_par_opt_vec ------------
    !--------------------- In case of power transformation, the structural parameters are back-transformed. ---------------------------------
    !------------------- At the end, d_S%struct_par_opt_vec, d_S%theta and d_S%sig are in the physical space. -------------------------------
    
    !Form d_S%struct_par_opt_vec
    d_S%struct_par_opt_vec = xmin !May include sigma
    
    !backtransform structural parameters
    k = 0
    if ((maxval(cv_S%struct_par_opt).eq.1)) then
        do i = 1, cv_PAR%p
            if (cv_S%struct_par_opt(i).eq.1) then
                if (cv_S%trans_theta(i).eq.1) then
                    do j = 1,cv_S%num_theta_type (i)
                        k = k + 1
                        d_S%struct_par_opt_vec(k) = ((d_S%struct_par_opt_vec(k) + cv_S%alpha_trans(i)) / (cv_S%alpha_trans(i)))**(cv_S%alpha_trans(i)) !Back-trans
                        d_S%theta(i,j) = d_S%struct_par_opt_vec(k)
                    enddo
                else
                    do j = 1,cv_S%num_theta_type (i)
                        k = k + 1
                        d_S%theta(i,j) = d_S%struct_par_opt_vec(k)
                    enddo
                endif
            endif
        enddo
    endif

    if (d_S%sig_opt.eq.1) then
        k=k+1
        if (d_S%trans_sig.eq.1) then
            d_S%struct_par_opt_vec(k) = ((d_S%struct_par_opt_vec(k) + d_S%alpha_trans_sig) / (d_S%alpha_trans_sig))**(d_S%alpha_trans_sig) !Back-trans
            d_S%sig = d_S%struct_par_opt_vec(k)
        else
            d_S%sig = d_S%struct_par_opt_vec(k)
        endif
    endif
    
    !*******************************************************************************************************************************************
    !Note: from here d_S%theta and d_S%sig contain the theta and sigma values optimized after the minimization of the structural parameter obj. function.
    !It's not true that the CURRENT Qss, Qsy, HQsy, Qyy and other variables, that depend on theta and sigma, are calculated with these optimized parameters.
    !This because in Nelder Mead the values that minimize the function can occur not during the last performed iteration.
    !*******************************************************************************************************************************************

    !Printing out the new structural parameters
    write (*,*) " "
    write (*,*) "The structural optimization found the following parameters"
    do p = 1,cv_PAR%p
        if (cv_S%var_type(Q0_All(p)%BetaAss) == 0) then
            imax = 1
        elseif (cv_S%var_type(Q0_All(p)%BetaAss) == 1 .or. cv_S%var_type(Q0_All(p)%BetaAss) == 2) then
            imax = 2
        elseif (cv_S%var_type(Q0_All(p)%BetaAss) == 3) then
            imax = cv_PAR%ndim
        endif
        do i = 1,imax
            write (*,334) p,i,d_S%theta(Q0_All(p)%BetaAss,i)
334         format(' Beta ass.',i2,', theta_0_',i1,': ',ES11.4)
        enddo
    enddo
    write (*,*) " "

    end subroutine PCGA_marginal_struct_param_optim

    !*************************************************
    !Evaluating the original structural parameters
    !*************************************************
    subroutine PCGA_beg_str_object_fun(d_XQR,cv_OBS,d_OBS,d_A,cv_S,d_PM,cv_PAR,cv_PM,d_S)

    implicit none
    type (kernel_XQR),    intent(in)     :: d_XQR
    type(cv_observ),      intent(in)     :: cv_OBS
    type(d_observ),       intent(in)     :: d_OBS
    type(d_algorithmic),  intent(inout)  :: d_A
    type(d_prior_mean),   intent(in)     :: d_PM
    type (d_struct),      intent(in)     :: d_S
    type(cv_param),       intent(in)     :: cv_PAR
    type(cv_struct),      intent(inout)  :: cv_S
    type(cv_prior_mean),  intent(in)     :: cv_PM

    integer                              :: errcode, i, curr_struct
    double precision, allocatable        :: z(:), pivot(:)
    double precision                     :: lndetGyy, ztiGyyz
    double precision, allocatable        :: UinvGyy(:,:) ! used as both U and InvGyy
    double precision, allocatable        :: TMPV(:)

    write (*,*) " "
    write (*,*) "Preparing for structural parameter optimization"

    allocate (pivot(cv_OBS%nobs))
    allocate (z(cv_OBS%nobs))
    allocate (UinvGyy(cv_OBS%nobs,cv_OBS%nobs))

    !----- intitialize variables
    lndetGyy  = 0.D0
    ztiGyyz   = 0.D0
    UinvGyy   = UNINIT_REAL   ! matrix

    ! At this point we need HXB, HXQbb, OMEGA only if we have prior information about beta.
    if (cv_PM%betas_flag.ne.0) then   !------> we have prior information about beta
        stop "Beta does not exist in PCGA bgaPEST. Therefore prior information does not estist"
    else !------> we don't have prior information about beta
        !Form the linearization residuals z
        z = d_OBS%obs + d_XQR%hs0Hs0 !d_XQR%hs0Hs0 has already a negative sign
    endif

    !*******************************************************************************
    !--- Calculate the determinant term of the str. pars objective function --------
    !----------------------------- 0.5ln(det(Gyy)) ---------------------------------
    !-- First perform LU decomposition on Gyy(=d_A%Qyy)
    UinvGyy = d_XQR%HQH !-nobs x nobs --- note that this is used as U in this context
    do i = 1, cv_OBS%nobs
        UinvGyy(i,i) = UinvGyy(i,i) + ((d_S%sig) *d_XQR%R0(i,i)) !adding R, diagonal.
    end do

    !LU decomposition
    call dgetrf(cv_OBS%nobs, cv_OBS%nobs, UinvGyy, cv_OBS%nobs, pivot, errcode) !MC: Error with pivot since it should be an integer. Error is handled by set the property Diagnostics > Check Routine Interfaces to No. This is however not not recommend: https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/282188.
    deallocate(pivot)
    do i = 1,cv_OBS%nobs
        lndetGyy = lndetGyy + dlog(abs(UinvGyy(i,i)))
    enddo
    lndetGyy = 0.5 * lndetGyy
    !*******************************************************************************

    !*******************************************************************************
    !------------ Calculate the misfit term of the objective function --------------
    !------------------------ 0.5(z' x invGyy x z) ---------------------------------
    !Calculate the inverse of Gyy(=d_A%Qyy)
    UinvGyy = d_XQR%HQH ! nobs x nobs, re-use UinvGyy, now as InvGyy
    do i = 1, cv_OBS%nobs
        UinvGyy(i,i) = UinvGyy(i,i) + ((d_S%sig) *d_XQR%R0(i,i)) !adding R, diagonal.
    end do

    call INVGM(cv_OBS%nobs,UinvGyy)
    !Form inv(Gyy)*z
    allocate(TMPV(cv_OBS%nobs))
    call DGEMV('n',cv_OBS%nobs, cv_OBS%nobs, 1.D0, UinvGyy, cv_OBS%nobs, &
        z, 1, 0.D0, TMPV,1) !On exit TMPV is invGyy*z
    !Multiply z' * TMPV and 0.5
    ! ddot ( n ,  x ,  incx ,  y ,  incy )
    ztiGyyz = 0.5 * ddot(cv_OBS%nobs, z, 1, TMPV, 1)
    if (allocated(TMPV)) deallocate(TMPV) !Deallocate TMPV no more necessary here
    !*******************************************************************************

    if (allocated(z))       deallocate(z)
    if (allocated(UinvGyy)) deallocate(UinvGyy)

    !******************************************************************************
    !----------------- OBJECTIVE FUNCTION FOR STRUCTURAL PARAMETERS ---------------
    cv_S%str_obj_fun = lndetGyy + ztiGyyz
    !******************************************************************************

    write (*,333) cv_S%str_obj_fun
333 format(' Objective function of the current structural parameters:',ES11.4)

    end subroutine PCGA_beg_str_object_fun

    !****************************************
    ! Nelder–Mead downhill simplex method
    !****************************************
    subroutine PCGA_nelmin_sp ( n, start, xmin, ynewlo, reqmin, step, konvge, kcount,icount, numres,ifault, &
        & d_XQR, Q0_all,cv_OBS, d_OBS, cv_A, d_A, d_PAR, cv_S, d_S, d_PM, cv_PAR,cv_PM,nQ0, errstruc, miostruc, cv_MIO, d_ANI, d_MIO, d_MOD)

    implicit none

    integer ( kind = 4 ) :: n
    real    ( kind = 8 ), parameter :: ccoeff = 0.5D+00
    real    ( kind = 8 ) :: del
    real    ( kind = 8 ) :: dn
    real    ( kind = 8 ) :: dnn
    real    ( kind = 8 ), parameter :: ecoeff = 2.0D+00
    real    ( kind = 8 ), parameter :: eps = 0.001D+00
    integer ( kind = 4 ) :: i
    integer ( kind = 4 ) :: icount
    integer ( kind = 4 ) :: ifault
    integer ( kind = 4 ) :: ihi
    integer ( kind = 4 ) :: ilo
    integer ( kind = 4 ) :: j
    integer ( kind = 4 ) :: jcount
    integer ( kind = 4 ) :: kcount
    integer ( kind = 4 ) :: konvge
    integer ( kind = 4 ) :: l
    integer ( kind = 4 ) :: nn
    integer ( kind = 4 ) :: numres
    real    ( kind = 8 ) :: p(n,n+1)
    real    ( kind = 8 ) :: p2star(n)
    real    ( kind = 8 ) :: pbar(n)
    real    ( kind = 8 ) :: pstar(n)
    real    ( kind = 8 ), parameter :: rcoeff = 1.0D+00
    real    ( kind = 8 ) :: reqmin
    real    ( kind = 8 ) :: rq
    real    ( kind = 8 ) :: start(n)
    real    ( kind = 8 ) :: step(n)
    real    ( kind = 8 ) :: x
    real    ( kind = 8 ) :: xmin(n)
    real    ( kind = 8 ) :: y(n+1)
    real    ( kind = 8 ) :: y2star
    real    ( kind = 8 ) :: ylo
    real    ( kind = 8 ) :: ynewlo
    real    ( kind = 8 ) :: ystar
    real    ( kind = 8 ) :: z

    integer                              :: nQ0
    type(kernel_XQR),     intent(inout)  :: d_XQR
    type(cv_observ),      intent(inout)  :: cv_OBS
    type(d_observ),       intent(inout)  :: d_OBS
    type(d_algorithmic),  intent(inout)  :: d_A
    type(cv_algorithmic), intent(inout)  :: cv_A
    type(d_prior_mean),   intent(in)     :: d_PM
    type(cv_param),       intent(inout)  :: cv_PAR
    type(d_param),        intent(inout)  :: d_PAR
    type(cv_struct),      intent(inout)  :: cv_S
    type(d_struct),       intent(inout)  :: d_S
    type (Q0_compr),      intent(inout), pointer :: Q0_All(:)
    type(cv_prior_mean),  intent(in)     :: cv_PM

    !Introduced for PCGA:
    type (err_failure_struc) , intent(inout)    :: errstruc
    type (mio_struc)         , intent(inout)    :: miostruc
    type (cv_minout)         , intent(in)       :: cv_MIO
    type (d_anisotropy)      , intent(in)       :: d_ANI
    type (d_minout)          , intent(in)       :: d_MIO
    type (d_comlin)          , intent(inout)    :: d_MOD

    !  Check the input parameters.
    if ( reqmin <= 0.0D+00 ) then
        ifault = 1
        return
    endif
    if ( n < 1 ) then
        ifault = 1
        return
    endif
    if ( konvge < 1 ) then
        ifault = 1
        return
    endif
    icount = 0
    numres = 0
    jcount = konvge
    dn = real ( n, kind = 8 )
    nn = n + 1
    dnn = real ( nn, kind = 8 )
    del = 1.0D+00
    rq = reqmin * dn
    !  Initial or restarted loop.
    do
        do i = 1, n
            p(i,nn) = start(i)
        enddo
        y(nn) = SPminN(start,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0, errstruc, miostruc, cv_MIO, d_ANI, d_MIO, d_MOD)
        icount = icount + 1
        do j = 1, n
            x = start(j)
            start(j) = start(j) + step(j) * del
            do i = 1, n
                p(i,j) = start(i)
            enddo
            y(j) = SPminN(start,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0, errstruc, miostruc, cv_MIO, d_ANI, d_MIO, d_MOD)
            icount = icount + 1
            start(j) = x
        enddo
        !  The simplex construction is complete.
        !  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
        !  the vertex of the simplex to be replaced.
        ylo = y(1)
        ilo = 1
        do i = 2, nn
            if ( y(i) < ylo ) then
                ylo = y(i)
                ilo = i
            endif
        enddo
        !  Inner loop.
        do
            if ( kcount <= icount ) then
                exit
            endif
            ynewlo = y(1)
            ihi = 1
            do i = 2, nn
                if ( ynewlo < y(i) ) then
                    ynewlo = y(i)
                    ihi = i
                endif
            enddo
            !  Calculate PBAR, the centroid of the simplex vertices
            !  excepting the vertex with Y value YNEWLO.
            do i = 1, n
                z = 0.0D+00
                do j = 1, nn
                    z = z + p(i,j)
                enddo
                z = z - p(i,ihi)
                pbar(i) = z / dn
            enddo
            !  Reflection through the centroid.
            do i = 1, n
                pstar(i) = pbar(i) + rcoeff * ( pbar(i) - p(i,ihi) )
            enddo
            ystar = SPminN(pstar,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0, errstruc, miostruc, cv_MIO, d_ANI, d_MIO, d_MOD)
            icount = icount + 1
            !  Successful reflection, so extension.
            if ( ystar < ylo ) then
                do i = 1, n
                    p2star(i) = pbar(i) + ecoeff * ( pstar(i) - pbar(i) )
                enddo
                y2star = SPminN(p2star,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0, errstruc, miostruc, cv_MIO, d_ANI, d_MIO, d_MOD)
                icount = icount + 1
                !  Check extension.
                if ( ystar < y2star ) then
                    do i = 1, n
                        p(i,ihi) = pstar(i)
                    enddo
                    y(ihi) = ystar
                    !  Retain extension or contraction.
                else
                    do i = 1, n
                        p(i,ihi) = p2star(i)
                    enddo
                    y(ihi) = y2star
                endif
                !  No extension.
            else
                l = 0
                do i = 1, nn
                    if ( ystar < y(i) ) then
                        l = l + 1
                    endif
                enddo
                if ( 1 < l ) then
                    do i = 1, n
                        p(i,ihi) = pstar(i)
                    enddo
                    y(ihi) = ystar
                    !  Contraction on the Y(IHI) side of the centroid.
                else if ( l == 0 ) then
                    do i = 1, n
                        p2star(i) = pbar(i) + ccoeff * ( p(i,ihi) - pbar(i) )
                    enddo
                    y2star = SPminN(p2star,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0, errstruc, miostruc, cv_MIO, d_ANI, d_MIO, d_MOD)
                    icount = icount + 1
                    !  Contract the whole simplex.
                    if ( y(ihi) < y2star ) then
                        do j = 1, nn
                            do i = 1, n
                                p(i,j) = ( p(i,j) + p(i,ilo) ) * 0.5D+00
                                xmin(i) = p(i,j)
                            enddo
                            y(j) = SPminN(xmin,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0, errstruc, miostruc, cv_MIO, d_ANI, d_MIO, d_MOD)
                            icount = icount + 1
                        enddo
                        ylo = y(1)
                        ilo = 1
                        do i = 2, nn
                            if ( y(i) < ylo ) then
                                ylo = y(i)
                                ilo = i
                            endif
                        enddo
                        cycle
                        !  Retain contraction.
                    else
                        do i = 1, n
                            p(i,ihi) = p2star(i)
                        enddo
                        y(ihi) = y2star
                    endif
                    !  Contraction on the reflection side of the centroid.
                else if ( l == 1 ) then
                    do i = 1, n
                        p2star(i) = pbar(i) + ccoeff * ( pstar(i) - pbar(i) )
                    enddo
                    y2star = SPminN(p2star,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0, errstruc, miostruc, cv_MIO, d_ANI, d_MIO, d_MOD)
                    icount = icount + 1
                    !  Retain reflection?
                    if ( y2star <= ystar ) then
                        do i = 1, n
                            p(i,ihi) = p2star(i)
                        enddo
                        y(ihi) = y2star
                    else
                        do i = 1, n
                            p(i,ihi) = pstar(i)
                        enddo
                        y(ihi) = ystar
                    endif
                endif
            endif
            !  Check if YLO improved.
            if ( y(ihi) < ylo ) then
                ylo = y(ihi)
                ilo = ihi
            endif
            jcount = jcount - 1
            if ( 0 < jcount ) then
                cycle
            endif
            !  Check to see if minimum reached.
            if ( icount <= kcount ) then
                jcount = konvge
                z = 0.0D+00
                do i = 1, nn
                    z = z + y(i)
                enddo
                x = z / dnn
                z = 0.0D+00
                do i = 1, nn
                    z = z + ( y(i) - x )**2
                enddo
                if ( z <= rq ) then
                    exit
                endif
            endif
        enddo
        !  Factorial tests to check that YNEWLO is a local minimum.
        do i = 1, n
            xmin(i) = p(i,ilo)
        enddo
        ynewlo = y(ilo)
        if ( kcount < icount ) then
            ifault = 2
            exit
        endif
        ifault = 0
        do i = 1, n
            del = step(i) * eps
            xmin(i) = xmin(i) + del
            z = SPminN(xmin,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0, errstruc, miostruc, cv_MIO, d_ANI, d_MIO, d_MOD)
            icount = icount + 1
            if ( z < ynewlo ) then
                ifault = 2
                exit
            endif
            xmin(i) = xmin(i) - del - del
            z = SPminN(xmin,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0, errstruc, miostruc, cv_MIO, d_ANI, d_MIO, d_MOD)
            icount = icount + 1
            if ( z < ynewlo ) then
                ifault = 2
                exit
            endif
            xmin(i) = xmin(i) + del
        enddo
        if ( ifault == 0 ) then
            exit
        endif
        !  Restart the procedure.
        do i = 1, n
            start(i) = xmin(i)
        enddo
        del = eps
        numres = numres + 1
    enddo
    return
    end subroutine PCGA_nelmin_sp

    !****************************************************************************************************************************
    !---------------- External function that contains the structural parameters objective function to minimize ------------------
    !****************************************************************************************************************************
    real (KIND = 8) function SPminN(str_par_opt_vec,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0, errstruc, miostruc, cv_MIO, d_ANI, d_MIO, d_MOD)

    implicit none
    integer                  , intent(in)               :: nQ0, n
    type (err_failure_struc) , intent(inout)            :: errstruc
    type (mio_struc)         , intent(inout)            :: miostruc
    type (cv_algorithmic)    , intent(inout)            :: cv_A
    type (cv_minout)         , intent(in)               :: cv_MIO
    type (cv_observ)         , intent(inout)            :: cv_OBS
    type (cv_param)          , intent(inout)            :: cv_PAR
    type (cv_prior_mean)     , intent(in)               :: cv_PM
    type (cv_struct)         , intent(inout)            :: cv_S
    type (d_anisotropy)      , intent(in)               :: d_ANI
    type (d_minout)          , intent(in)               :: d_MIO
    type (d_comlin)          , intent(inout)            :: d_MOD
    type (d_observ)          , intent(inout)            :: d_OBS
    type (d_param)           , intent(inout)            :: d_PAR
    type (d_prior_mean)      , intent(in)               :: d_PM
    type (d_struct)          , intent(inout)            :: d_S
    type (kernel_XQR)        , intent(inout)            :: d_XQR
    type (Q0_compr)          , intent(inout), pointer   :: Q0_All(:)
    double precision         , intent(inout)            :: str_par_opt_vec (n)!This must be the pars vector to be optimized for. Must be the first argument in PCGA_SP_min (NelMead requires this)

    double precision         , allocatable              :: z(:), pivot(:), dtheta(:), TMPV(:), UinvGyy(:,:) !used as both U and InvGyy
    double precision                                    :: lndetGyy, ztiGyyz, dthQttdth
    integer                                             :: errcode, i, j, k, imax

    allocate (z(cv_OBS%nobs))
    allocate (pivot(cv_OBS%nobs))
    allocate (UinvGyy(cv_OBS%nobs,cv_OBS%nobs))

    !----- intitialize variables
    errcode   = UNINIT_INT
    lndetGyy  = 0.D0
    ztiGyyz   = 0.D0
    dthQttdth = 0.D0
    UinvGyy   = UNINIT_REAL   ! matrix

    !******************************************************************************************************
    !First we need to reform d_S%theta and d_S%sig overwriting the elements that must be optimized for and
    !leaving unchanged the others. These elements are into str_par_opt_vec (passed by NelMead).
    !d_S%theta and d_S%sig are used during the matrix operations.(Qss Qsy HQsy Qyy)
    !*****************************************************************************************************
    !In case of power transformation str_par_opt_vec may contains values in the estimation space. The value
    !will be back-transformed and assigned to d_S%struct_par_opt_vec.
    !*****************************************************************************************************
    d_S%struct_par_opt_vec = str_par_opt_vec
    k = 0
    if ((maxval(cv_S%struct_par_opt).eq.1)) then
        do i = 1, cv_PAR%p
            if (cv_S%struct_par_opt(i).eq.1) then
                if (cv_S%trans_theta(i).eq.1) then
                    do j = 1,cv_S%num_theta_type (i)
                        k = k + 1
                        d_S%struct_par_opt_vec(k) = ((str_par_opt_vec(k) + cv_S%alpha_trans(i)) / (cv_S%alpha_trans(i)))**(cv_S%alpha_trans(i))
                        d_S%theta(i,j) = d_S%struct_par_opt_vec(k)
                    enddo
                else
                    do j = 1,cv_S%num_theta_type (i)
                        k = k + 1
                        d_S%theta(i,j) = d_S%struct_par_opt_vec(k)
                    enddo
                endif
            endif
        enddo
    endif

    if (d_S%sig_opt.eq.1) then
        k=k+1
        if (d_S%trans_sig.eq.1) then
            d_S%struct_par_opt_vec(k) = ((str_par_opt_vec(k) + d_S%alpha_trans_sig) / (d_S%alpha_trans_sig))**(d_S%alpha_trans_sig)
            d_S%sig = d_S%struct_par_opt_vec(k)
        else
            d_S%sig = d_S%struct_par_opt_vec(k)
        endif
    endif
    !****************************************************************************************************************************
    !****************************************************************************************************************************

    !Print out parameters which are tested
    write (*,*) " "
    write (*,*) 'Testing:'
    do j = 1,cv_PAR%p
        if (cv_S%var_type(Q0_All(j)%BetaAss) == 0) then
            imax = 1
        elseif (cv_S%var_type(Q0_All(j)%BetaAss) == 1 .or. cv_S%var_type(Q0_All(j)%BetaAss) == 2) then
            imax = 2
        elseif (cv_S%var_type(Q0_All(j)%BetaAss) == 3) then
            imax = cv_PAR%ndim
        endif
        do i = 1,imax
            write (*,334) j,i,d_S%theta(Q0_All(j)%BetaAss,i)
334         format(' Beta ass.',i2,', theta_0_',i1,': ',ES11.4)
        enddo
        write (*,*) " "
    enddo

    !****************************************************************************************************************************
    ! At this point d_S%theta, d_S%sig and d_S%struct_par_opt_vec are ready to be used in the calculation of the obj func.
    !---------------------------------- d_S%struct_par_opt_vec is in the physical space. ----------------------------------------
    !****************************************************************************************************************************
    !We need to recalculate Qyy = H*Qss*H^T + R
    !First constructing a new low rank approximation of Qss if theta opimization is required.
    if ((maxval(cv_S%struct_par_opt).eq.1)) then
        call Low_rank_app(cv_A, cv_PAR, cv_S, d_ANI, d_PAR, d_S, d_XQR, Q0_All)
    endif

    !Now we compute the new Qyy:
    !Transforming the parameters to parameters-estimation space
    if (maxval(d_PM%Partrans).ge.1) then  !If yes, the parameters transformation is required
        call PCGA_sen_par_trans(d_PAR%pars, cv_PAR, d_PAR, d_PM) !Converting sensitivity and parameters in the estimation space
    end if

    !Forward computation
    call matrix_free_approach_computations(miostruc,errstruc,cv_A,cv_MIO,cv_PAR,d_PAR,cv_OBS,d_MIO,d_OBS,d_MOD,d_XQR,d_PM,2)
    
    !Transforming the parameters to the physical space
    if (maxval(d_PM%Partrans).ge.1) then  !If yes, we need to back-transform the parameters in the physical space
        call PCGA_par_back_trans(d_PAR%pars, cv_PAR, d_PAR, d_PM) !Converting sensitivity and parameters in the estimation space
    end if

    !Compute BC
    allocate(d_XQR%BC(cv_OBS%nobs,cv_A%K))
    !dgemm ( transa='n' ,  transb='n' ,  m=cv_OBS%nobs ,  n=cv_A%K ,  k=cv_A%K ,  alpha=1.d0 ,  a=d_XQR%B ,  lda=cv_OBS%nobs ,  b=d_XQR%tC ,  ldb=cv_A%K ,  beta=0.d0 ,  c=d_XQR%BC ,  ldc=cv_OBS%nobs )
    call dgemm ( 'n' ,  'n' ,  cv_OBS%nobs ,  cv_A%K ,  cv_A%K ,  1.d0 ,  d_XQR%B ,  cv_OBS%nobs ,  d_XQR%tC ,  cv_A%K ,  0.d0 ,  d_XQR%BC ,  cv_OBS%nobs )

    !Compute HQH^T = BC*B^T
    if (associated(d_XQR%HQH)) deallocate(d_XQR%HQH)
    allocate(d_XQR%HQH(cv_OBS%nobs,cv_OBS%nobs))
    !dgemm ( transa='n' ,  transb='t' ,  m=cv_OBS%nobs ,  n=cv_OBS%nobs ,  k=cv_A%K ,  alpha=1.d0 ,  a=d_XQR%BC ,  lda=cv_OBS%nobs ,  b=B ,  ldb=cv_OBS%nobs ,  beta=0.d0 ,  c=d_XQR%HQH ,  ldc=cv_OBS%nobs )
    call dgemm ( 'n' ,  't' ,  cv_OBS%nobs ,  cv_OBS%nobs ,  cv_A%K ,  1.d0 ,  d_XQR%BC ,  cv_OBS%nobs ,  d_XQR%B ,  cv_OBS%nobs ,  0.d0 ,  d_XQR%HQH ,  cv_OBS%nobs )

    deallocate(d_XQR%B)
    deallocate(d_XQR%BC)

    ! At this point we need HXB, HXQbb, OMEGA only if we have prior information about beta.
    if (cv_PM%betas_flag.ne.0) then   !------> we have prior information about beta
        stop "Beta does not exist in PCGA bgaPEST. Therefore prior information about this parameter does not exist"
    else !------> we don't have prior information about beta
        !Form the linearization residuals z
        z = d_OBS%obs + d_XQR%hs0Hs0
    endif
    !***********************************************************************************************************

    !*******************************************************************************
    !--- Calculate the determinant term of the str. pars objective function --------
    !----------------------------- 0.5ln(det(Gyy)) ---------------------------------
    !-- First perform LU decomposition on Gyy
    UinvGyy = d_XQR%HQH !-nobs x nobs --- note that this is used as U in this context
    do i = 1, cv_OBS%nobs
        UinvGyy(i,i) = UinvGyy(i,i) + ((d_S%sig) *d_XQR%R0(i,i)) !adding R, diagonal.
    end do

    call dgetrf(cv_OBS%nobs, cv_OBS%nobs, UinvGyy, cv_OBS%nobs, pivot, errcode)
    if (allocated(pivot))       deallocate(pivot)
    do i = 1,cv_OBS%nobs
        lndetGyy = lndetGyy + dlog(abs(UinvGyy(i,i)))
    enddo
    lndetGyy = 0.5 * lndetGyy
    !*******************************************************************************

    !*******************************************************************************
    !------------ Calculate the misfit term of the objective function --------------
    !------------------------ 0.5(z' x invGyy x z) ---------------------------------
    !Calculate the inverse of Gyy
    UinvGyy = d_XQR%HQH   ! nobs x nobs, re-use UinvGyy, now as InvGyy
    do i = 1, cv_OBS%nobs
        UinvGyy(i,i) = UinvGyy(i,i) + ((d_S%sig) *d_XQR%R0(i,i)) !adding R, diagonal.
    end do

    call INVGM(cv_OBS%nobs,UinvGyy)
    !Form inv(Gyy)*z
    allocate(TMPV(cv_OBS%nobs))
    call DGEMV('n',cv_OBS%nobs, cv_OBS%nobs, 1.D0, UinvGyy, cv_OBS%nobs, &
        z, 1, 0.D0, TMPV,1) !On exit TMPV is invGyy*z
    !Multiply z' * TMPV and 0.5
    ztiGyyz = 0.5 * ddot(cv_OBS%nobs, z, 1, TMPV, 1)
    !call DGEMV('t',cv_OBS%nobs, 1, 5.0D-1, z, cv_OBS%nobs, TMPV, 1, 0.D0, ztiGyyz,1)
    if (allocated(TMPV)) deallocate(TMPV) !Deallocate TMPV no more necessary here
    !*******************************************************************************

    !***************************************************************************************************
    !------------------- Calculate the prior theta/sig term of the objective function ------------------
    !--------------------------------- 0.5(dtheta x invQtheta x dtheta) --------------------------------
    !--> note: if theta covariance form is set as 0 (meaning no prior covariance on theta provided) and
    !--- sig_p_var is 0 too, assume totally unknown and do not consider dthQttdth in calculations.
    !-- This "if" statement is not strictly necessary (because with the previous assumptions dthQttdth
    !-- will be zero anyway), but we avoid useless computations.
    if  ((cv_A%theta_cov_form.ne.0).or.(d_S%sig_p_var.ne.0.)) then !"if" not strictly necessary, see above
        allocate (dtheta(cv_S%num_theta_opt))
        dtheta    = UNINIT_REAL   ! matrix
        !Form dtheta
        dtheta = d_S%struct_par_opt_vec - d_S%struct_par_opt_vec_0
        !Form invQtt*dtheta
        allocate(TMPV(cv_S%num_theta_opt))
        call DGEMV('n',cv_S%num_theta_opt, cv_S%num_theta_opt, 1.D0, d_S%invQtheta, &
            cv_S%num_theta_opt, dtheta, 1, 0.D0, TMPV,1) !On exit TMPV is invQtt*dtheta
        !Multiply dtheta' * TMPV and 0.5
        dthQttdth = 0.5 * ddot(cv_S%num_theta_opt, dtheta, 1, TMPV, 1)
        !call DGEMV('t',cv_S%num_theta_opt, 1, 5.0D-1, dtheta, cv_S%num_theta_opt, &
        !    TMPV, 1, 0.D0, dthQttdth,1)
        if (allocated(TMPV))   deallocate(TMPV) !Deallocate TMPV no more necessary here
        if (allocated(dtheta)) deallocate(dtheta)
    endif
    !**************************************************************************************************

    if (allocated(z))       deallocate(z)
    if (allocated(UinvGyy)) deallocate(UinvGyy)

    !******************************************************************************
    !----------------- OBJECTIVE FUNCTION FOR STRUCTURAL PARAMETERS ---------------
    SPminN = lndetGyy + ztiGyyz + dthQttdth
    !******************************************************************************
    
    !Printing the objective function of the tested structural parameters
    write(*,333) SPminN
333 format(' Current test of structural parameters gave the objective function:',ES11.4)
    return

    end function SPminN
    end module PCGA_struct_param_optimization