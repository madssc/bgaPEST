    module PCGA_normalize_X

    use bayes_pest_control
    use utilities
    use error_message
    use jupiter_input_data_support

    contains

    !***********************************************************************
    ! This function finds the normalized version of X and saves it in U_X
    !***********************************************************************
    subroutine Normalize_X(d_XQR, cv_PAR)

    implicit none
    type(kernel_XQR), intent(inout) :: d_XQR
    type(cv_param),   intent(in)    :: cv_PAR
    double precision, allocatable   :: work_X(:), Sigma_X(:)
    integer                         :: work_X_size, info_X
    character (len=ERRORWIDTH)      :: retmsg

    !****************
    ! Normalize X
    !****************
    nullify(d_XQR%U_X)
    allocate(d_XQR%U_X(cv_PAR%npar,cv_PAR%p))

    !Number of columns in X
    if (cv_PAR%p == 1) then
        d_XQR%U_X(:,:) = d_XQR%X(:,:) / sqrt(dble(cv_PAR%npar)) !dble = double()

        !If there are several beta associations then X will have several columns.
        !In the future X could also be a function of positions.
    else if (cv_PAR%p > 1) then
        allocate(Sigma_X(cv_PAR%p))
        
        d_XQR%U_X(:,:) = d_XQR%X(:,:)
        
        !Using SVD to normalize X. First using a test to find the propper parameters.
        !dgesvd ( jobu='O' ,  jobvt='N' ,  m=cv_PAR%npar ,  n=cv_PAR%p ,  a=d_XQR%U_X ,  lda=cv_PAR%npar ,  s=Sigma_X ,  u=unused 0 ,  ldu=cv_PAR%npar ,  vt=unused 0,  ldvt=1 ,  work=bestwork ,  lwork=-1 ,  info=info_X )
        allocate(work_X(1))
        call dgesvd('O', 'N', cv_PAR%npar, cv_PAR%p, d_XQR%U_X, cv_PAR%npar, Sigma_X, 0, cv_PAR%npar, 0, cv_PAR%p, work_X, -1, info_X)
        work_X_size = int(work_X(1))
        deallocate(work_X)
        allocate(work_X(work_X_size))
        call dgesvd('O', 'N', cv_PAR%npar, cv_PAR%p, d_XQR%U_X, cv_PAR%npar, Sigma_X, 0.d0, cv_PAR%npar, 0.d0, cv_PAR%p, work_X, work_X_size, info_X)
        deallocate(work_X)
        deallocate(Sigma_X)
    else
        write(retmsg,12)
12      format('Warning: cv_PAR%p is less than 1 i.e. X is has zero width.')
        call utl_writmess(6,retmsg)
        stop
    endif

    if (allocated(d_XQR%X))     deallocate(d_XQR%X)
    
    end subroutine Normalize_X
    end module PCGA_normalize_X