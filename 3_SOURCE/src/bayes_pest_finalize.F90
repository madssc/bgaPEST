    module bayes_pest_finalize

    contains
    subroutine bpd_finalize(estruc, mstruc, Q0_All, d_PM, d_S, cv_A, cv_S, cv_PAR, d_PAR, &
        &    cv_OBS, d_OBS, d_MOD, d_MIO, d_XQR)
    use bayes_pest_control
    use error_message
    use model_input_output
    implicit none
    !--  Main Data Arrays for OBS and PARS

    type(err_failure_struc), intent(inout)  :: estruc
    type(mio_struc),         intent(inout)  :: mstruc
    type (Q0_compr),         pointer        :: Q0_All(:)
    type (d_prior_mean)  :: d_PM
    type (d_algorithmic) :: d_A
    type (kernel_XQR)    :: d_XQR
    type (d_struct)      :: d_S
    type (cv_algorithmic):: cv_A
    type (cv_struct)     :: cv_S
    type (cv_param)      :: cv_PAR
    type (d_param)       :: d_PAR
    type (cv_observ)     :: cv_OBS
    type (d_observ)      :: d_OBS
    type (d_comlin)      :: d_MOD
    type (d_minout)      :: d_MIO
    integer              :: ierr, mio_err, p, status

    !-- Clean up the mio
    mio_err=mio_finalise(estruc,mstruc)
    if (mio_err /= 0) write(*,*) "Mio structure not deallocate correct"

    !-- deallocate Q0_all
    do p = 1,cv_PAR%p
        if (associated(Q0_All(p)%Q0_C)) deallocate(Q0_All(p)%Q0_C,stat=ierr)
        if (cv_A%derivative_test /= 1) then
            if (associated(Q0_All(p)%tU_Z)) deallocate(Q0_All(p)%tU_Z,stat=ierr)
            if (associated(Q0_All(p)%tC))   deallocate(Q0_All(p)%tC,stat=ierr)
        endif
    enddo
    if (associated(Q0_All)) deallocate(Q0_All,stat=ierr)

    !-- deallocate algorithmic data array (d_A)
    if (associated(d_A%H))           deallocate(d_A%H,stat=ierr)
    if (associated(d_A%HQHt))        deallocate(d_A%HQHt,stat=ierr)
    if (associated(d_A%Hsold))       deallocate(d_A%Hsold,stat=ierr)
    if (associated(d_A%Qsy))         deallocate(d_A%Qsy,stat=ierr)
    if (associated(d_A%Qyy))         deallocate(d_A%Qyy,stat=ierr)

    !-- deallocate prior mean information (d_PM)
    if (associated(d_PM%beta_0))     deallocate(d_PM%beta_0,stat=ierr)
    if (associated(d_PM%Qbb))        deallocate(d_PM%Qbb,stat=ierr)
    if (associated(d_PM%InvQbb))     deallocate(d_PM%InvQbb,stat=ierr)
    if (associated(d_PM%InvQbbB0))   deallocate(d_PM%InvQbbB0,stat=ierr)
    if (associated(d_PM%alpha))      deallocate(d_PM%alpha,stat=ierr)
    if (associated(d_PM%Partrans))   deallocate(d_PM%Partrans,stat=ierr)

    !-- deallocate structural parameter data (d_XQR)
    if (associated(d_XQR%X))         deallocate(d_XQR%X,stat=ierr)
    if (associated(d_XQR%Q0))        deallocate(d_XQR%Q0,stat=ierr)
    if (associated(d_XQR%R0))        deallocate(d_XQR%R0,stat=ierr)

    !PCGA:
    if (allocated(d_XQR%BC))      deallocate(d_XQR%BC,stat=ierr)
    if (allocated(d_XQR%X))       deallocate(d_XQR%X,stat=ierr)
    if (allocated(d_XQR%U_X))     deallocate(d_XQR%U_X,stat=ierr)
    if (allocated(d_XQR%Ap))      deallocate(d_XQR%Ap,stat=ierr)
    if (allocated(d_XQR%HQ))      deallocate(d_XQR%HQ,stat=ierr)
    if (allocated(d_XQR%LA_p))    deallocate(d_XQR%LA_p,stat=ierr)
    if (allocated(d_XQR%MA_p))    deallocate(d_XQR%MA_p,stat=ierr)
    if (associated(d_XQR%Ap))     deallocate(d_XQR%Ap,stat=ierr)
    if (associated(d_XQR%HQ))     deallocate(d_XQR%HQ,stat=ierr)
    if (associated(d_XQR%LA_p))   deallocate(d_XQR%LA_p,stat=ierr)
    if (associated(d_XQR%MA_p))   deallocate(d_XQR%MA_p,stat=ierr)
    !if (associated(d_XQR%mu))     deallocate(d_XQR%MA_p,stat=ierr)
    if (associated(d_XQR%tU_Z))   deallocate(d_XQR%tU_Z,stat=ierr)
    if (associated(d_XQR%tC))     deallocate(d_XQR%tC,stat=ierr)
    if (associated(d_XQR%U_X))    deallocate(d_XQR%U_X,stat=ierr)
    if (associated(d_XQR%invtC))  deallocate(d_XQR%invtC,stat=ierr)
    if (associated(d_XQR%hs0Hs0)) deallocate(d_XQR%hs0Hs0,stat=ierr)
    if (associated(d_XQR%HQH))    deallocate(d_XQR%HQH,stat=ierr)
    if (associated(d_XQR%LApy))   deallocate(d_XQR%LApy,stat=ierr)
    if (cv_A%PCGA == 1) then
        status = vsldeletestream ( cv_A%stream_Gaussian )
    endif
    
    !-- deallocate structural parameter data (d_S)
    if (associated(d_S%theta_0))     deallocate(d_S%theta_0,stat=ierr)
    if (associated(d_S%theta_cov))   deallocate(d_S%theta_cov,stat=ierr)
    if (associated(d_S%invQtheta))   deallocate(d_S%invQtheta,stat=ierr)
    if (associated(d_S%struct_par_opt_vec_0)) deallocate(d_S%struct_par_opt_vec_0,stat=ierr)
    if (associated(d_S%struct_par_opt_vec))   deallocate(d_S%struct_par_opt_vec,stat=ierr)
    if (associated(d_S%theta))       deallocate(d_S%theta,stat=ierr)

    !-- deallocate parameter control values (cv_PAR)
    if (associated(cv_PAR%grp_name))   deallocate(cv_PAR%grp_name,stat=ierr)
    if (associated(cv_PAR%grp_type))   deallocate(cv_PAR%grp_type,stat=ierr)
    if (associated(cv_PAR%derinc))     deallocate(cv_PAR%derinc,stat=ierr)
    if (associated(cv_PAR%K))          deallocate(cv_PAR%K,stat=ierr)
    if (associated(cv_PAR%K_s))        deallocate(cv_PAR%K_s,stat=ierr)

    !-- deallocate parameter structure (d_PAR)
    if (associated(d_PAR%group))              deallocate(d_PAR%group,stat=ierr)
    if (associated(d_PAR%pars))               deallocate(d_PAR%pars,stat=ierr)
    if (associated(d_PAR%pars_perturbed))     deallocate(d_PAR%pars_perturbed,stat=ierr) !Mads
    if (associated(d_PAR%parnme))             deallocate(d_PAR%parnme,stat=ierr)
    if (associated(d_PAR%pars_old))           deallocate(d_PAR%pars_old,stat=ierr)
    if (associated(d_PAR%pars_lns))           deallocate(d_PAR%pars_lns,stat=ierr)
    if (associated(d_PAR%lox))                deallocate(d_PAR%lox,stat=ierr)
    if (associated(d_PAR%SenMethod))          deallocate(d_PAR%SenMethod,stat=ierr)
    if (associated(d_PAR%BetaAssoc))          deallocate(d_PAR%BetaAssoc,stat=ierr)
    if (associated(d_PAR%Group_type))         deallocate(d_PAR%Group_type,stat=ierr)

    !-- deallocate observation structure (d_OBS)
    if (associated(d_OBS%group))     deallocate(d_OBS%group,stat=ierr)
    if (associated(d_OBS%obs))       deallocate(d_OBS%obs,stat=ierr)
    if (associated(d_OBS%obsnme))    deallocate(d_OBS%obsnme,stat=ierr)
    if (associated(d_OBS%h))         deallocate(d_OBS%h,stat=ierr)
    if (associated(d_OBS%weight))    deallocate(d_OBS%weight,stat=ierr)
    if (associated(d_OBS%TemplateFileNo)) deallocate(d_OBS%TemplateFileNo,stat=ierr)

    !-- deallocate observation (cv_OBS)
    if (associated(cv_OBS%grp_name)) deallocate(cv_OBS%grp_name,stat=ierr)

    !-- deallocate model i/o structure (d_MIO)
    if (associated(d_MIO%tpl))       deallocate(d_MIO%tpl,stat=ierr)
    if (associated(d_MIO%infle))     deallocate(d_MIO%infle,stat=ierr)
    if (associated(d_MIO%ins))       deallocate(d_MIO%ins,stat=ierr)
    if (associated(d_MIO%outfle))    deallocate(d_MIO%outfle,stat=ierr)
    if (associated(d_MIO%pargroup))  deallocate(d_MIO%pargroup,stat=ierr)
    if (associated(d_MIO%derinc_PCGA)) deallocate(d_MIO%derinc_PCGA,stat=ierr)

    !-- deallocate (cv_S)
    if (associated(cv_S%prior_cov_mode)) deallocate(cv_S%prior_cov_mode,stat=ierr)
    if (associated(cv_S%var_type))       deallocate(cv_S%var_type,stat=ierr)
    if (associated(cv_S%struct_par_opt)) deallocate(cv_S%struct_par_opt,stat=ierr)
    if (associated(cv_S%num_theta_type)) deallocate(cv_S%num_theta_type,stat=ierr)
    if (associated(cv_S%trans_theta))    deallocate(cv_S%trans_theta,stat=ierr)
    if (associated(cv_S%alpha_trans))    deallocate(cv_S%alpha_trans,stat=ierr)

    end subroutine bpd_finalize

    end module bayes_pest_finalize