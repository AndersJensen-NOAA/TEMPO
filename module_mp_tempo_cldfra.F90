! Prognostic cloud fraction module for TEMPO Microphysics
!=================================================================================================================
module module_mp_tempo_cldfra
    use module_mp_tempo_params, only : lvap0, lsub, cp2, R, Rv, R1, critical_rh, cf_low, D0r, am_r, am_i
    use module_mp_tempo_utils, only : rslf, rsif

#if defined(mpas)
    use mpas_kind_types, only: wp => RKIND, sp => R4KIND, dp => R8KIND
#elif defined(standalone)
    use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
#else
    use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
#define ccpp_default 1
#endif
#if defined(ccpp_default) && defined(MPI)
    use mpi_f08
#endif

    implicit none
    private
    public:: tempo_cldfra_driver

  contains
!=================================================================================================================    
  subroutine tempo_cldfra_driver(i,j,kts,kte,dt,temp,pres,rho,w,qv,qa,qabl,qc,qcbl,qi,qt,qcexp,qr,nr, &
       qiexp,qs,ncexp,niexp)
!=================================================================================================================

    integer, intent(in) :: i, j, kts, kte
    real, intent(in) :: dt
    real, intent(in) :: pres(:), rho(:), w(:), qt(:), qabl(:), qcbl(:)
    real, intent(inout) :: qa(:), qv(:), qc(:), temp(:), qcexp(:), qr(:), nr(:), qi(:), qiexp(:), qs(:), &
         ncexp(:), niexp(:)

    integer :: k
    real :: tempc
    real, parameter :: p0 = 100000.
!    real, parameter :: ls_w_limit = 1.
    real, parameter :: ls_w_limit = 0.5
    real, parameter :: grav = 9.8
    real :: al, bs, sd, qc_calc
    real :: term1, term2, term3, gterm, eros_term
    real :: cf_check, qrten, qsten, lami, xtmp, lam0
    logical, dimension(kts:kte) :: create_sgs_clouds, erode_sgs_clouds, evolve_sgs_clouds
    real, dimension(kts:kte) :: qvs, qvi, U, U00, lvap, ocp, dqsdT, omega, lsub2, lfus2
    real, dimension(kts:kte) :: qaten_create, qaten_evolve, qaten_erode
    real, dimension(kts:kte) :: qcten_create, qcten_evolve, qcten_erode
    real, dimension(kts:kte) :: thten_create, thten_evolve, thten_erode
    real, dimension(kts:kte) :: dcond_ls

    real, dimension(kts:kte) :: qaten, qcten, qiten, thten, thten_l, thten_i, qaten_l

    do k = kts, kte

       if ((qc(k) <= R1) .and. (qt(k) <= R1)) then
          qc(k) = qc(k) + qcbl(k)
          qa(k) = max(qa(k), qabl(k))
       endif

       create_sgs_clouds(k) = .false.
       erode_sgs_clouds(k) = .false.
       evolve_sgs_clouds(k) = .false.
       dcond_ls(k) = 0.
       qaten_create(k) = 0.
       qcten_create(k) = 0.
       thten_create(k) = 0.
       qaten_evolve(k) = 0.
       qcten_evolve(k) = 0.
       thten_evolve(k) = 0.
       qaten_erode(k) = 0.
       qcten_erode(k) = 0.
       thten_erode(k) = 0.
       qaten(k) = 0.
       qcten(k) = 0.
       qiten(k) = 0.
       thten(k) = 0.

       ! Environment
       tempc = temp(k) - 273.15
       qvs(k) = rslf(pres(k), temp(k))
       if (tempc <= 0.0) then
          qvi(k) = rsif(pres(k), temp(k))
       else
          qvi(k) = qvs(k)
       endif
       U(k) = qv(k) / qvs(k)
       U00(k) = critical_rh

       lvap(k) = lvap0 + (2106.0 - 4218.0)*tempc
       lsub2(k) = lsub
       lfus2(k) = lsub2(k) - lvap(k)
       ocp(k) = 1.0 / (cp2*(1.+0.887*qv(k)))
       omega(k) = -grav * w(k) * rho(k)
       dqsdT(k) = lvap(k) * qvs(k) / (Rv*temp(k)**2)

       ! Limit cloud fraction to range of cf_low (in module_mp_tempo_params) to 100%
       if ((qc(k) <= R1) .and. (qr(k) <= R1) .and. (qi(k) <= R1)) then
          qa(k) = 0.0
       else
          qa(k) = min(qa(k), 1.0)
          qa(k) = max(qa(k), cf_low)
       endif

       ! Create liquid clouds if no cloud water or ice present, otherwise, evolve SGS clouds
       if((qa(k) < 1.0) .and. (w(k) < ls_w_limit) .and. (U(k) < 1.0) .and. (U(k) > U00(k))) then
          if (qc(k) <= R1) then
             if ((qt(k) <= R1) .and. (w(k) > 0.) .and. (tempc > -30.)) then
                create_sgs_clouds(k) = .true.
             endif
          else
             evolve_sgs_clouds(k) = .true.
          endif
       endif

       ! Erosion
       if(qa(k) < 1.0 .and. (qc(k) > R1) .and. (qt(k) <= R1) .and. U(k) < 1.) then
          erode_sgs_clouds(k) = .true.
       endif

       ! Varibles used in cloud fraction scheme
       al = 1. / (1. + dqsdT(k)*lvap(k)*ocp(k))
       bs = al * (1.-U00(k)) * qvs(k)
       sd = al*(qvs(k)-qv(k))
       qc_calc = al*(qv(k)+qc(k)-qvs(k))

       ! Either create or evolve clouds
       ! Wilson and Gregory 2003, Eqs. 27, 28
       if (create_sgs_clouds(k)) then
          qaten_create(k) = 0.5/bs*(bs+qc_calc)/dt
          qcten_create(k) = qaten_create(k)*0.5*(bs+qc_calc)
          thten_create(k) = ((p0/pres(k))**(R/cp2))*lvap(k)*ocp(k)*qcten_create(k)
          ! Create at least 5% cloud fraction
          if ((qa(k) + qaten_create(k)*dt) < 0.05 .or. (qc(k) + qcten_create(k)*dt) <= R1) then
             qcten_create(k) = 0.
             qaten_create(k) = 0.
             thten_create(k) = 0.
          endif
       elseif (evolve_sgs_clouds(k)) then
          ! Contains large-scale forcing and longwave cooling
          ! Closure terms from Wilson and Gregory 2003, Eq. 22
          dcond_ls(k) = -al * dqsdT(k) * (omega(k)/rho(k)*ocp(k)) ! + lwrad(k)*(pres(k)/p0)**(R/cp2))
          if ((abs(sd) > R1) .and. (dcond_ls(k) > 0.) .and. (qc(k) > 0.)) then
             term1 = ((1.-qa(k))**2*qa(k)**2/qc(k)) + (qa(k)**2*(1.-qa(k))**2/sd)
             term2 = (1.-qa(k))**2 + qa(k)**2
             gterm = 0.5*term1/term2
             qaten_evolve(k) = gterm*dcond_ls(k)
             qcten_evolve(k) = qa(k)*dcond_ls(k)
             thten_evolve(k) = ((p0/pres(k))**(R/cp2))*lvap(k)*ocp(k)*qcten_evolve(k)
          else
             qcten_evolve(k) = 0.
             qaten_evolve(k) = 0.
             thten_evolve(k) = 0.
          endif
       endif

       ! Erosion
       ! Wilson et al. 2008, Eqs. A11, A12
       if (erode_sgs_clouds(k)) then
          if ((abs(sd) > R1) .and. (qc(k) > 0.)) then
             term1 = ((1.-qa(k))**2*qa(k)**2/qc(k)) + (qa(k)**2*(1.-qa(k))**2/sd)
          else
             term1 = 0.
          endif
          term2 = (1.-qa(k))**2 + qa(k)**2
          gterm = 0.5*term1/term2
          term3 = -3.1*qc_calc/(al*qvs(k))
          eros_term = -2.25e-5 * exp(term3)
          qcten_erode(k) = min(0., ((qc(k) - qc_calc*qa(k))*eros_term))
          qaten_erode(k) = min(0., (-gterm*qc_calc*eros_term))
          thten_erode(k) = ((p0/pres(k))**(R/cp2))*lvap(k)*ocp(k)*qcten_erode(k)
       endif

       ! Sum tendencies
       qaten(k) = qaten_create(k) + qaten_evolve(k) + qaten_erode(k)
       qcten(k) = qcten_create(k) + qcten_evolve(k) + qcten_erode(k)
       thten(k) = thten_create(k) + thten_evolve(k) + thten_erode(k)

       ! If eroding cloud fraction to less than 5%, or cloud water
       ! to less than R1, then completely remove
!        if ((((qa(k) + qaten(k)*dt) < 0.05) .and. (qaten(k) < 0.)) .or. &
!             ((qc(k) + qcten(k)*dt) <= R1)) then
!            qaten(k) = -qa(k)/dt
!            qcten(k) = -qc(k)/dt
!            thten(k) = ((p0/pres(k))**(R/cp2))*lvap(k)*ocp(k)*qcten(k)
!         endif

!        ! Update variables
!        qc(k) = qc(k) + qcten(k)*dt
!        qv(k) = qv(k) - qcten(k)*dt
!        qa(k) = qa(k) + qaten(k)*dt
!        temp(k) = temp(k) + ((thten(k)) / ((p0/pres(k))**(R/cp2)))*dt

!        if ((qc(k) > R1)) then
!           if ((qa(k) > 0.98) .and. (qcexp(k) > R1)) then
!              qcexp(k) = qcexp(k) + qc(k)
! !!!             ncexp(k) = ncexp(k) + 10.e6 ! 10 cm^-3
!              qc(k) = 0.
! ! could still be rain             qa(k) = 0.
!           else
!              qrten = min((1350.*qc(k)**2.47*(10.)**-1.79), qc(k)/dt)
!              qr(k) = qr(k) + qrten*dt
!              nr(k) = nr(k) + qrten*dt / (am_r*(4.*D0r)**3.0)
!              qc(k) = qc(k) - qrten*dt
!           endif
!        endif

!        ! Final check on cloud fraction
!        if ((qc(k) <= R1) .and. (qr(k) <= R1)) then
!           qa(k) = 0.0
!        else
!           qa(k) = min(qa(k), 1.0)
!           qa(k) = max(qa(k), cf_low)
!        endif

!!    enddo

    ! Ice
!!    do k = kts, kte
       create_sgs_clouds(k) = .false.
       erode_sgs_clouds(k) = .false.
       evolve_sgs_clouds(k) = .false.
       dcond_ls(k) = 0.
       qaten_create(k) = 0.
       qcten_create(k) = 0.
       thten_create(k) = 0.
       qaten_evolve(k) = 0.
       qcten_evolve(k) = 0.
       thten_evolve(k) = 0.
       qaten_erode(k) = 0.
       qcten_erode(k) = 0.
       thten_erode(k) = 0.
!       qaten(k) = 0.
!       qcten(k) = 0.
!       thten(k) = 0.

       ! Environment
!       tempc = temp(k) - 273.15
!       qvs(k) = rslf(pres(k), temp(k))
!       if (tempc <= 0.0) then
!          qvi(k) = rsif(pres(k), temp(k))
!       else
!          qvi(k) = qvs(k)
!       endif
       U(k) = qv(k) / qvi(k)
!       U00(k) = critical_rh

       ! ocp(k) = 1.0 / (cp2*(1.+0.887*qv(k)))
       ! omega(k) = -grav * w(k) * rho(k)
       dqsdT(k) = lsub2(k) * qvi(k) / (Rv*temp(k)**2)

       ! Limit cloud fraction to range of cf_low (in module_mp_tempo_params) to 100%
 !      if ((qc(k) <= R1) .and. (qr(k) <= R1)) then
 !         qa(k) = 0.0
 !      else
 !         qa(k) = min(qa(k), 1.0)
 !         qa(k) = max(qa(k), cf_low)
 !      endif

       ! Create liquid clouds if no cloud water or ice present, otherwise, evolve SGS clouds
       if((qa(k) < 1.0) .and. (w(k) < ls_w_limit) .and. (U(k) < 1.0) .and. (U(k) > U00(k))) then
          if (qi(k) <= R1) then
             if ((w(k) > 0.) .and. (tempc <= -30.) .and. (qt(k) <= R1) .and. (qc(k) <= R1)) then
                create_sgs_clouds(k) = .true.
             endif
          else
             if (tempc < 0.) evolve_sgs_clouds(k) = .true.
          endif
       endif

       ! Erosion
       if(qa(k) < 1.0 .and. (qi(k) > R1) .and. (qt(k) <= R1) .and. U(k) < 1.) then
          erode_sgs_clouds(k) = .true.
       endif

       ! Varibles used in cloud fraction scheme
       al = 1. / (1. + dqsdT(k)*lsub2(k)*ocp(k))
       bs = al * (1.-U00(k)) * qvi(k)
       sd = al*(qvi(k)-qv(k))
       qc_calc = al*(qv(k)+qi(k)-qvi(k))

       ! Either create or evolve clouds
       ! Wilson and Gregory 2003, Eqs. 27, 28
       if (create_sgs_clouds(k)) then
          qaten_create(k) = 0.5/bs*(bs+qc_calc)/dt
          qcten_create(k) = qaten_create(k)*0.5*(bs+qc_calc)
          thten_create(k) = ((p0/pres(k))**(R/cp2))*lsub2(k)*ocp(k)*qcten_create(k)
          ! Create at least 5% cloud fraction
          if ((qa(k) + qaten_create(k)*dt) < 0.05 .or. (qi(k) + qcten_create(k)*dt) <= R1) then
             qcten_create(k) = 0.
             qaten_create(k) = 0.
             thten_create(k) = 0.
          endif
       elseif (evolve_sgs_clouds(k)) then
          ! Contains large-scale forcing and longwave cooling
          ! Closure terms from Wilson and Gregory 2003, Eq. 22
          dcond_ls(k) = -al * dqsdT(k) * (omega(k)/rho(k)*ocp(k)) ! + lwrad(k)*(pres(k)/p0)**(R/cp2))
          if ((abs(sd) > R1) .and. (dcond_ls(k) > 0.) .and. (qi(k) > 0.)) then
             term1 = ((1.-qa(k))**2*qa(k)**2/qi(k)) + (qa(k)**2*(1.-qa(k))**2/sd)
             term2 = (1.-qa(k))**2 + qa(k)**2
             gterm = 0.5*term1/term2
             qaten_evolve(k) = gterm*dcond_ls(k)
             qcten_evolve(k) = qa(k)*dcond_ls(k)
             thten_evolve(k) = ((p0/pres(k))**(R/cp2))*lsub2(k)*ocp(k)*qcten_evolve(k)
          else
             qcten_evolve(k) = 0.
             qaten_evolve(k) = 0.
             thten_evolve(k) = 0.
          endif
       endif

       ! Erosion
       ! Wilson et al. 2008, Eqs. A11, A12
       if (erode_sgs_clouds(k)) then
          if ((abs(sd) > R1) .and. (qi(k) > 0.)) then
             term1 = ((1.-qa(k))**2*qa(k)**2/qi(k)) + (qa(k)**2*(1.-qa(k))**2/sd)
          else
             term1 = 0.
          endif
          term2 = (1.-qa(k))**2 + qa(k)**2
          gterm = 0.5*term1/term2
          term3 = -3.1*qc_calc/(al*qvi(k))
          eros_term = -2.25e-5 * exp(term3)
          qcten_erode(k) = min(0., ((qi(k) - qc_calc*qa(k))*eros_term))
          qaten_erode(k) = min(0., (-gterm*qc_calc*eros_term))
          thten_erode(k) = ((p0/pres(k))**(R/cp2))*lsub2(k)*ocp(k)*qcten_erode(k)
       endif

       ! Sum tendencies
       qaten(k) = qaten(k) + qaten_create(k) + qaten_evolve(k) + qaten_erode(k)
       qiten(k) = qcten_create(k) + qcten_evolve(k) + qcten_erode(k)
       thten(k) = thten(k) + thten_create(k) + thten_evolve(k) + thten_erode(k)

       ! If eroding cloud fraction to less than 5%, or cloud water
       ! to less than R1, then completely remove
       if (((qa(k) + qaten(k)*dt) < 0.05) .and. (qaten(k) < 0.)) then
          qaten(k) = -qa(k)/dt
          qcten(k) = -qc(k)/dt
          qiten(k) = -qi(k)/dt
          thten(k) = ((p0/pres(k))**(R/cp2))*lvap(k)*ocp(k)*qcten(k) + &
               ((p0/pres(k))**(R/cp2))*lsub2(k)*ocp(k)*qiten(k)
       else
          if ((qc(k) + qcten(k)*dt) <= R1) then
             qcten(k) = -qc(k)/dt
             thten(k) = ((p0/pres(k))**(R/cp2))*lvap(k)*ocp(k)*qcten(k)
          endif
          if ((qi(k) + qiten(k)*dt) <= R1) then
             qiten(k) = -qi(k)/dt
             thten(k) = thten(k) + ((p0/pres(k))**(R/cp2))*lsub2(k)*ocp(k)*qiten(k)
          endif
       endif

       ! Update variables
       qc(k) = qc(k) + qcten(k)*dt
       qi(k) = qi(k) + qiten(k)*dt
       qv(k) = qv(k) - qcten(k)*dt - qiten(k)*dt
       qa(k) = qa(k) + qaten(k)*dt
       temp(k) = temp(k) + ((thten(k)) / ((p0/pres(k))**(R/cp2)))*dt

       ! Conversions
       if ((qc(k) > R1)) then
          if ((qa(k) > 0.98) .and. (qcexp(k) > R1)) then
             qcten(k) = qc(k)/dt
             qcexp(k) = qcexp(k) + qcten(k)*dt
             ! Keep low (or zero out)
             ncexp(k) = ncexp(k) + qcten(k)*dt / (am_r*(50.e-6)**3.0)
             qc(k) = qc(k) - qcten(k)*dt
          else
             qrten = min((1350.*qc(k)**2.47*(10.)**(-1.79)), qc(k)/dt)
             qr(k) = qr(k) + qrten*dt
             nr(k) = nr(k) + qrten*dt / (am_r*(4.*D0r)**3.0)
             qc(k) = qc(k) - qrten*dt
          endif
          if (tempc < -30.) then
             qcten(k) = qc(k)/dt
             qi(k) = qi(k) + qcten(k)*dt
             qc(k) = qc(k) - qcten(k)*dt
             thten(k) = +((p0/pres(k))**(R/cp2))*lfus2(k)*ocp(k)*qcten(k)/dt
             temp(k) = temp(k) + ((thten(k)) / ((p0/pres(k))**(R/cp2)))*dt
          endif
       endif

       if ((qi(k) > R1)) then
          if ((qa(k) > 0.98) .and. (qiexp(k) > R1)) then
             qiten(k) = qi(k)/dt
             qiexp(k) = qiexp(k) + qiten(k)*dt
             niexp(k) = niexp(k) + qiten(k)*dt / (am_i*(150.e-6)**3.0)
             qi(k) = qi(k) - qiten(k)*dt
          else
             lami = (am_i*6.*10.e3/rho(k)/qi(k))**(0.33333)
             lam0 = 5000.
             qsten = 0.
             if (lami < lam0) then
                xtmp = lami/lam0
                qsten = qi(k)/dt*min((1-xtmp**3),1.)
             endif
             qs(k) = qs(k) + qsten*dt
             qi(k) = qi(k) - qsten*dt
          endif
          if (tempc > 0.) then
             qiten(k) = qi(k)/dt
             qc(k) = qc(k) + qiten(k)*dt
             qi(k) = qi(k) - qiten(k)*dt
             thten(k) = -((p0/pres(k))**(R/cp2))*lfus2(k)*ocp(k)*qiten(k)/dt
             temp(k) = temp(k) + ((thten(k)) / ((p0/pres(k))**(R/cp2)))*dt
          endif
       endif

       ! Final check on cloud fraction
       if ((qc(k) <= R1) .and. (qr(k) <= R1) .and. (qi(k) <= R1)) then
          qa(k) = 0.0
       else
          qa(k) = min(qa(k), 1.0)
          qa(k) = max(qa(k), cf_low)
       endif

    enddo

  end subroutine tempo_cldfra_driver

!=================================================================================================================    
end module module_mp_tempo_cldfra
!=================================================================================================================
