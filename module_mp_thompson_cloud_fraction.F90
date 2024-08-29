module module_mp_thompson_cloud_fraction

  implicit none
  
contains

  subroutine subgrid_scale_liquid_clouds(kts, kte, no_micro, l_qc, dt, temp, rho, pres, lvap, ocp, qvs, w1d, lw1d, &
       qa, qv, rc, nc, rr, prw_sgi, prw_sgb, prw_sge)

    integer, intent(in) :: kts, kte
    real, intent(in) :: rho(:), pres(:), lvap(:), ocp(:), qvs(:), w1d(:), lw1d(:)
    real, intent(in) :: dt
    logical, intent(inout) :: no_micro, l_qc(:)
    real, intent(inout) :: temp(:), qa(:), qv(:), rc(:), nc(:), rr(:)
    double precision, intent(inout) :: prw_sgi(:), prw_sgb(:), prw_sge(:)
    
    real :: rvgas = 461.50
    real, dimension(kts:kte) :: dcond_ls, qa_tend_ls, cooling_ls, dqsdT, U, U00, omega
    real, dimension(kts:kte) :: dcond_ls_eros, qa_tend_ls_eros, rc_tmp
    logical, dimension(kts:kte) :: do_sgs_clouds, create_clouds, erode_sgs_clouds
    real :: al, bs, sd, qc_calc, term1, term2, gterm, eros_term, term3, ratio
    integer :: k
    
    do k = kts, kte
       
!       if ((k >= 6) .and. (qa(k) > 0.)) then
!           write(*,*) 'BEGIN', k, qa(k), rc(k), prw_sge(k), prw_sgb(k), prw_sgi(k), U(k)
!       endif
       
       dcond_ls(k) = 0.
       qa_tend_ls(k) = 0.
       dcond_ls_eros(k) = 0.
       qa_tend_ls_eros(k) = 0.
       cooling_ls(k) = 0.
       do_sgs_clouds(k) = .false.
       erode_sgs_clouds(k) = .false.
       create_clouds(k) = .false.
       dqsdT(k) = lvap(k) * qvs(k) / (rvgas*temp(k)**2)
       U(k) = qv(k) / qvs(k)
       U00(k) = 0.9
       omega(k) = -9.8 * w1d(k) * rho(k)
!!       write(*,*) k, lw1d(k)
       
       if (rc(k)*qa(k)/rho(k) <= 1.e-12) then
          qa(k) = 0.0
       else
          qa(k) = min(qa(k), 1.0)
          qa(k) = max(qa(k), 0.01)
       endif
       
       if(qa(k) < 1.0 .and. rc(k)*qa(k)/rho(k) > 1.e-12 .and. (U(k) < 1.) .and. &
            w1d(k) < 0.1 .and. (U(k) > (U00(k) + 0.01))) then
          do_sgs_clouds(k) = .true.
       endif
       
       !if(qa(k) < 1.0 .and. rc(k)*qa(k)/rho(k) > 1.e-12 .and. (U(k) < 1.) .and. &
       !     (U(k) > (U00(k) + 0.01))) then
       !   do_sgs_clouds(k) = .true.
       !endif

       if(qa(k) < 1.0 .and. rc(k)*qa(k)/rho(k) > 1.e-12 .and. U(k) < 1.) then
          erode_sgs_clouds(k) = .true.
       endif

       ! ADD CHECK ON RI TOO (DONT MAKE SGS IF RI EXISTS)
       !if(qa(k) < 1. .and. rc(k)*qa(k)/rho(k) <= 1.e-12 .and. (U(k) < 1.) .and. &
       !     w1d(k) < 0.1 .and. rr(k)/rho(k) <= 1.e-12 .and. (U(k) > (U00(k) + 0.01))) then
       !   create_clouds(k) = .true.
       !endif

       ! THIS WORKS BETTER
       if(qa(k) < 1. .and. rc(k)*qa(k)/rho(k) <= 1.e-12 .and. (U(k) < 1.) .and. &
            w1d(k) < 0.1 .and. (U(k) > (U00(k) + 0.01))) then
          create_clouds(k) = .true.
       endif
       
       !if(qa(k) < 1. .and. rc(k)*qa(k)/rho(k) <= 1.e-12 .and. (U(k) < 1.) .and. &
       !     (U(k) > (U00(k) + 0.01))) then
       !   create_clouds(k) = .true.
       !endif

    enddo

    do k = kts, kte
       al = 1. / (1. + dqsdT(k)*lvap(k)*ocp(k))
       bs = al * (1.-U00(k)) * qvs(k)
       sd = al*(qvs(k)-qv(k))
       qc_calc = al*(qv(k)+(rc(k)*qa(k)/rho(k))-qvs(k))

       if (create_clouds(k)) then
          qa_tend_ls(k) = 0.5/bs*(bs+qc_calc)/dt
          dcond_ls(k) = qa_tend_ls(k)*0.5*(bs+qc_calc)

          ! Update the state
          if ((rc(k)*qa(k)/rho(k)+dcond_ls(k)*dt) > 1.e-12) then
             no_micro = .false.
             l_qc(k) = .true.
             qa(k) = qa(k) + qa_tend_ls(k)*dt
             qa(k) = min(qa(k), 1.0)
             qa(k) = max(qa(k), 0.01)
             ! prw_sgi(k) = dcond_ls(k)*rho(k)*qa(k)
             prw_sgi(k) = dcond_ls(k)/qa(k)
             
!             if ((qa(k) + qa_tend_ls(k)*dt) < 0.05) then
!                prw_sgi(k) = 0.
!             else
!                no_micro = .false.
!                l_qc(k) = .true.
!                qa(k) = qa(k) + qa_tend_ls(k)*dt
!                qa(k) = min(qa(k), 1.0)
!                qa(k) = max(qa(k), 0.01)
!                prw_sgi(k) = dcond_ls(k)*rho(k)*qa(k)
!             endif

!             write(*,*) 'CREATE CLOUDS', k, qa(k), rc(k), prw_sgi(k), U(k), qa_tend_ls(k), dcond_ls(k)

             
             ! Create at least 5% or some reasonable amount of cloud water (0.001 g m-3)
             ! Create at least 5%
             if (qa(k) < 0.05) then
             !! if ((qa(k) < 0.05) .or. ((rc(k)/rho(k)+dcond_ls(k)*dt) < 1.e-6)) then
                prw_sgi(k) = 0.0
                no_micro = .true.
                l_qc(k) = .false.
                qa(k) = 0.
             endif
             
          endif
       elseif (do_sgs_clouds(k)) then
          dcond_ls(k) = -al * dqsdT(k) * (omega(k)/rho(k)*ocp(k) + lw1d(k)*(pres(k)/100000.)**0.286 ) * qa(k)
          if ((abs(sd) > 1.e-12) .and. (dcond_ls(k) > 0.)) then
             term1 = (1.-qa(k))**2*qa(k)/(rc(k)/rho(k)) + qa(k)**2*(1-qa(k))**2/sd
             term2 = (1.-qa(k))**2 + qa(k)**2
             gterm = 0.5*term1/term2
             qa_tend_ls(k) = gterm/qa(k)*dcond_ls(k)
          else
             dcond_ls(k) = 0.
             qa_tend_ls(k) = 0.
          endif

          ! Update the state
!          if ((rc(k)*qa(k)/rho(k)+dcond_ls(k)*dt) > 1.e-12 .and. (qa(k)+qa_tend_ls(k)*dt > 0.05)) then
          if ((rc(k)*qa(k)/rho(k)+dcond_ls(k)*dt) > 1.e-12) then
             no_micro = .false.
             l_qc(k) = .true.

!             if ((qa(k) > 0.1) .and. (qa(k) < 0.9)) then
!                write(*,*) 'BUILD CLOUDS 1', k, qa(k), rc(k), prw_sgb(k), U(k), qa_tend_ls(k), dcond_ls(k)
!             endif
             qa(k) = qa(k) + qa_tend_ls(k)*dt
!             prw_sgb(k) = dcond_ls(k)*rho(k)*qa(k)
             qa(k) = min(qa(k), 1.0)
             qa(k) = max(qa(k), 0.01)
             prw_sgb(k) = dcond_ls(k)/qa(k)
             
             

!             if ((qa(k) > 0.1) .and. (qa(k) < 0.9)) then
!                write(*,*) 'BUILD CLOUDS 2', k, qa(k), rc(k), prw_sgb(k), U(k), qa_tend_ls(k), dcond_ls(k)
!             endif
             
          endif
       endif

       if (erode_sgs_clouds(k)) then
          no_micro = .false.
          l_qc(k) = .true.
          if ((abs(sd) > 1.e-12)) then
             term1 = (1.-qa(k))**2*qa(k)/(rc(k)/rho(k)) + qa(k)**2*(1-qa(k))**2/sd
          else
             term1 = 0.
          endif
          
          term2 = (1.-qa(k))**2 + qa(k)**2
          gterm = 0.5*term1/term2
          
          ! Erosion
          term3 = -3.1*qc_calc/(al*qvs(k))
          eros_term = -2.25e-5 * exp(term3)
          
          dcond_ls(k) = min(0., (rc(k)*qa(k)/rho(k) - qc_calc*qa(k))*eros_term)
          qa_tend_ls(k) = min(0., -gterm*qc_calc*eros_term)

!          if ((qa(k) > 0.1) .and. (qa(k) < 0.9)) then
!             write(*,*) 'ERODE CLOUDS 1', k, qa(k), rc(k), prw_sge(k), U(k), qa_tend_ls(k), dcond_ls(k)
!          endif
          
          ! Update the state
          qa(k) = qa(k) + qa_tend_ls(k)*dt
!!          prw_sge(k) = max(qa(k)*dcond_ls(k), -rc(k)/rho(k)/dt)

          qa(k) = min(qa(k), 1.0)
          qa(k) = max(qa(k), 0.01)
          
          prw_sge(k) = max(dcond_ls(k)/qa(k), -rc(k)/rho(k)/dt)

!          if ((qa(k) > 0.1) .and. (qa(k) < 0.9)) then
!             write(*,*) 'ERODE CLOUDS 2', k, qa(k), rc(k), prw_sge(k), U(k), qa_tend_ls(k), dcond_ls(k)
!          endif
          !          if (qa_tend_ls(k) < -1.) then

          
          ! Erode sgs cloud field at some point
!          if (qa(k) < 0.05 .and. rc(k)*qa(k) < 1.e-9 .and. (U(k) <= (U00(k) + 0.01))) then
!!!          if (qa(k) < 0.025 .and. (U(k) <= (U00(k) + 0.01))) then
!             prw_sge(k) = -rc(k)/rho(k)/dt
!             prw_sgi(k) = 0.
!             prw_sgb(k) = 0.
!             qa(k) = 0.0
!          endif
          
          
!          if (qa(k) < 0.05 .and. rc(k) < 1.e-9 .and. (U(k) <= (U00(k) + 0.01))) then
!             prw_sge(k) = -rc(k)/rho(k)/dt
!             prw_sgi(k) = 0.
!             prw_sgb(k) = 0.
!             qa(k) = 0.0
!          else
!             qa(k) = min(qa(k), 1.0)
!             qa(k) = max(qa(k), 0.01)
!          endif
          
       endif
       
       ! Test clean up
       ! if (qa(k) < 0.05 .and. rc(k) <= 1.e-12 .and. (U(k) <= (U00(k) + 0.01))) then
       !    qa(k) = 0.0
       ! endif
       
!       if ((k >= 6) .and. (qa(k) > 0.)) then
!          write(*,*) 'END', k, qa(k), rc(k), prw_sge(k), prw_sgb(k), prw_sgi(k), U(k)
!       endif
    enddo
       

    
  end subroutine subgrid_scale_liquid_clouds

  
  subroutine init_sgs_clouds(temp, pres, qa, qv, qvs, ls_or_lv, ocp, &
       dcond_ls, qa_tend, cooling_terms, dqsdT, gamm, U, U00, do_sgs_clouds)

    real, dimension(:), intent(in) :: temp, pres, qa, qv, qvs, ls_or_lv, ocp
    real, dimension(:), intent(inout) :: dcond_ls, qa_tend, cooling_terms, dqsdT, gamm, U
    logical, intent(inout) :: do_sgs_clouds(:)
    
    integer :: kte, k
    real :: rvgas = 461.50
    real :: U00, alpha1, al, bs
    
    kte = size(dcond_ls)
    do k = 1, kte
       dcond_ls(k) = 0.
       qa_tend(k) = 0.
       cooling_terms(k) = 0.
       do_sgs_clouds(k) = .false.
       dqsdT(k) = ls_or_lv(k) * qvs(k) / (rvgas*temp(k)**2)
       gamm(k) = dqsdT(k) * ls_or_lv(k) * ocp(k)
       U(k) = qv(k) / qvs(k)
       U00 = 0.85

       if((qa(k) > 0.) .and. (qa(k) < 1.0)) then
          do_sgs_clouds(k) = .true.
       endif

       if((qa(k) == 0.) .and. (U(k) > (U00 + 0.01))) then
          do_sgs_clouds(k) = .true.
       endif
       
    enddo
    
  end subroutine init_sgs_clouds

  
  subroutine subgrid_scale_ice_clouds(kts, kte, dt, rho, pres, lvap, lsub, ocp, qvs, qvsi, w1d, &
       no_micro_grid, qa, qv, rc, nc, ri, ni, temp)

    ! Input / Output
    integer, intent(in) :: kts, kte
    real, intent(in) :: rho(:), pres(:), ocp(:), qvs(:), qvsi(:), w1d(:)
    real, intent(in) :: dt, lvap(:), lsub
    logical, intent(in) :: no_micro_grid(:)
    real, intent(inout) :: qa(:), qv(:), rc(:), nc(:), ri(:), ni(:), temp(:)

    ! Constants
    real :: con_eps = 0.608
    real :: con_rd = 287.
    real :: grav = 9.8
    real :: cp_air = 1004.
    real :: hlv = 2.5e6

    ! Local variables
    integer :: k
    real :: dcond_ls_l(kts:kte), qa_tend_l(kts:kte), qc_tend(kts:kte)
    real :: dcond_ls_i(kts:kte), qa_tend_i(kts:kte)

    real :: dqsdT_l(kts:kte), dqsdT_i(kts:kte), lsub_k(kts:kte), ls_or_lv(kts:kte)
    real :: gamm_l(kts:kte), U_l(kts:kte), omega(kts:kte)
    real :: gamm_i(kts:kte), U_i(kts:kte), U00_l
    real :: dqs_ls_l(kts:kte), dqs_ls_i(kts:kte), tmp1(kts:kte), tmp2(kts:kte), A_dt(kts:kte), B_dt(kts:kte)
    real :: qa0(kts:kte), qaeq(kts:kte), qa1(kts:kte), qabar(kts:kte), da_ls(kts:kte)
    real :: cooling_terms_l(kts:kte), cooling_terms_i(kts:kte)
    real :: eros(kts:kte)
    real :: term1, term2, term3, term4, sink, ratio_l, ratio_i, pressure_threshold
    logical :: do_sgs_clouds_l(kts:kte), do_sgs_clouds_i(kts:kte)
    real :: al, bs, gterm, sd, qcdt, qc_calc
    
    do k = kts, kte
       omega(k) = -grav * w1d(k) * rho(k)
       lsub_k(k) = lsub
    enddo
    
    ! Liquid
    ! call init_sgs_clouds(temp, pres, qa, qv, qvs, lvap, ocp, &
    !      dcond_ls_l, qa_tend_l, cooling_terms_l, dqsdT_l, gamm_l, U_l, U00_l, ls_or_lv, do_sgs_clouds_l)

    ! do k = kts, kte
    !    if(do_sgs_clouds_l(k)) then
          
    !       al = 1. / (1. + dqsdT_l(k)*lvap(k)*ocp(k))
    !       bs = al * (1.-U00_l) * qvs(k)
    !       sd = al*(qvs(k)-qv(k))
          
    !       if (qa(k) == 0.) then
    !          ! Init clouds
    !          qc_calc = al*(qv(k)-qvs(k))
    !          qa_tend_l(k) = 0.5*bs*(bs+qc_calc)/dt
    !          qc_tend(k) = qa_tend_l(k)*0.5*(bs+qc_calc)/dt
    !       else
    !          !Evolve clouds
    !          term1 = 0.0
    !          if (rc(k) > 1.e-12) then
    !             term1 = (1.-qa(k))**2*qa(k)/rc(k) + qa(k)**2*(1-qa(k))**2/sd
    !          endif
    !          term2 = (1.-qa(k))**2 + qa(k)**2
    !          gterm = 0.5*term1/term2

    !          qa_tend_l(k) = -al * dqsdT_l(k) * (omega(k)/rho(k)*ocp(k)) * gterm
    !          qc_tend(k) = -al * dqsdT_l(k) * (omega(k)/rho(k)*ocp(k)) * qa(k)

    !          ! Update the state
    !          qa(k) = qa(k) + qa_tend_l(k)*dt
    !          rc(k) = rc(k) + qc_tend(k)*dt
    !          nc(k) = nc(k) + ((qc_tend(k)*dt) / (3.14159*1000./6.*((10.e-6)**3.0)))
       
    !          qv(k) = qv(k) - qc_tend(k)*dt
    !          temp(k) = temp(k) + lvap(k)*ocp(k)*qc_tend(k)*dt

    !       endif

    !       ! Erosion
          
          
    !    endif
    ! enddo
          
  end subroutine subgrid_scale_ice_clouds
  
end module module_mp_thompson_cloud_fraction
