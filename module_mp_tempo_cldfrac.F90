! Cloud fraction module for TEMPO Microphysics
!=================================================================================================================
module module_mp_tempo_cldfrac
    use module_mp_tempo_params, only : lvap0, lsub, cp2, R, Rv, R1, critical_rh, cf_low
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
    public:: tempo_cldfrac_driver
    
  contains
!=================================================================================================================    
  subroutine tempo_cldfrac_driver(i,j,kts,kte,dt,temp,pres,rho,w,lwrad,qa,qv,qc,qi,nc,ni)
!=================================================================================================================
    
    integer, intent(in) :: i, j, kts, kte
    real, intent(in) :: dt
    real, intent(in) :: pres(:), rho(:), w(:), lwrad(:)
    real, intent(inout) :: qa(:), qv(:), qc(:), qi(:), nc(:), ni(:), temp(:)
    
    integer :: k
    real :: tempc
    real, parameter :: p0 = 100000.
    real, parameter :: ls_w_limit = 1.
    real, parameter :: grav = 9.8
    real :: al, ali, bs, bsi, sd, sdi, qc_calc, qi_calc
    real :: term1, term2, term3, gterm, eros_term
    real :: cf_check
    real, dimension(kts:kte) :: qaten, qcten, thten, qiten, thten_l, thten_i, qaten_l, qaten_i
    real, dimension(kts:kte) :: qvs, qvi, U, Ui, U00, U00i, lvap, ocp, dqsdT, dqidT, omega
    real, dimension(kts:kte) :: qcten_create, qcten_evolve, qcten_erode
    real, dimension(kts:kte) :: qaten_create, qaten_evolve, qaten_erode
    real, dimension(kts:kte) :: thten_create, thten_evolve, thten_erode
    real, dimension(kts:kte) :: qiten_create, qiten_evolve, qiten_erode
    real, dimension(kts:kte) :: qaiten_create, qaiten_evolve, qaiten_erode
    real, dimension(kts:kte) :: thiten_create, thiten_evolve, thiten_erode
    real, dimension(kts:kte) :: cooling_ls, dcond_ls, lw1d, dcondi_ls
    logical, dimension(kts:kte) :: create_sgs_clouds, erode_sgs_clouds, evolve_sgs_clouds
    logical, dimension(kts:kte) :: create_sgs_ice_clouds, erode_sgs_ice_clouds, evolve_sgs_ice_clouds

    do k = kts, kte
       create_sgs_clouds(k) = .false.
       erode_sgs_clouds(k) = .false.
       evolve_sgs_clouds(k) = .false.
       create_sgs_ice_clouds(k) = .false.
       erode_sgs_ice_clouds(k) = .false.
       evolve_sgs_ice_clouds(k) = .false.
       qaten(k) = 0.
       qaten_l(k) = 0.
       qaten_i(k) = 0.
       qcten(k) = 0.
       qiten(k) = 0.
       thten(k) = 0.
       thten_l(k) = 0.
       thten_i(k) = 0.
       qaten_create(k) = 0.
       qcten_create(k) = 0.
       thten_create(k) = 0.
       qaten_evolve(k) = 0.
       qcten_evolve(k) = 0.
       thten_evolve(k) = 0.
       qaten_erode(k) = 0.
       qcten_erode(k) = 0.
       thten_erode(k) = 0.
       qaiten_create(k) = 0.
       qiten_create(k) = 0.
       thiten_create(k) = 0.
       qaiten_evolve(k) = 0.
       qiten_evolve(k) = 0.
       thiten_evolve(k) = 0.
       qaiten_erode(k) = 0.
       qiten_erode(k) = 0.
       thiten_erode(k) = 0.
       dcond_ls(k) = 0.
       dcondi_ls(k) = 0.
       lw1d(k) = 0.
       cooling_ls(k) = 0.
       tempc = temp(k) - 273.15
       qvs(k) = rslf(pres(k), temp(k))
!       qvi(k) = rsif(pres(k), temp(k))
       if (tempc <= 0.0) then
          qvi(k) = rsif(pres(k), temp(k))
       else
          qvi(k) = qvs(k)
       endif
       U(k) = qv(k) / qvs(k)
       Ui(k) = qv(k) / qvi(k)
       U00(k) = critical_rh
       U00i(k) = critical_rh

       lvap(k) = lvap0 + (2106.0 - 4218.0)*tempc
       ocp(k) = 1.0 / (cp2*(1.+0.887*qv(k)))
       omega(k) = -grav * w(k) * rho(k)
       dqsdT(k) = lvap(k) * qvs(k) / (Rv*temp(k)**2)
       dqidT(k) = lsub * qvi(k) / (Rv*temp(k)**2)

       ! Limit cloud fraction to range of cf_low to 100%
       ! if cloud water present
       if ((qc(k) <= R1) .and. (qi(k) <= R1)) then
          qa(k) = 0.0
       else
          qa(k) = min(qa(k), 1.0)
          qa(k) = max(qa(k), cf_low)
       endif

       ! Create SGS liquid clouds if no cloud water present
       ! Otherwise, evolve SGS clouds from forcing and LW cooling
       if((qa(k) < 1.0) .and. (w(k) < ls_w_limit) .and. (U(k) < 1.0) .and. (U(k) > U00(k))) then
          if (qc(k) <= R1) then
             if ((w(k) > 0.) .and. (qi(k) <= R1) .and. (tempc > -20.)) then
                create_sgs_clouds(k) = .true.
             endif
          else
             ! qc present so evolve qc sgs clouds
             evolve_sgs_clouds(k) = .true.
          endif
       endif

       if((qa(k) < 1.0) .and. (w(k) < ls_w_limit) .and. (Ui(k) < 1.0) .and. (Ui(k) > U00i(k))) then
          if (qi(k) <= R1) then
             if ((w(k) > 0.) .and. (qc(k) <= R1) .and. (tempc <= -20.)) then
!!                create_sgs_ice_clouds(k) = .true.
             endif
          else
             ! qc present so evolve qc sgs clouds
!!             evolve_sgs_ice_clouds(k) = .true.
          endif
       endif

       ! Erode SGS clouds
       if(qa(k) < 1.0 .and. (qc(k) > R1) .and. U(k) < 1.) then
          erode_sgs_clouds(k) = .true.
       endif

!       if(qa(k) < 1.0 .and. (qi(k) > R1) .and. Ui(k) < 1.) then
!          erode_sgs_ice_clouds(k) = .true.
!       endif

       ! Varibles used in cloud fraction scheme
       al = 1. / (1. + dqsdT(k)*lvap(k)*ocp(k))
       ali = 1. / (1. + dqidT(k)*lsub*ocp(k))
       bs = al * (1.-U00(k)) * qvs(k)
       bsi = ali * (1.-U00i(k)) * qvi(k)
       sd = al*(qvs(k)-qv(k))
       sdi = ali*(qvi(k)-qv(k))
       qc_calc = al*(qv(k)+qc(k)-qvs(k))
       qi_calc = ali*(qv(k)+qi(k)-qvi(k))

       ! Either create or evolve clouds
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
          ! Contains large-scale forcing and LW cooling (not currently used)
          dcond_ls(k) = -al * dqsdT(k) * (omega(k)/rho(k)*ocp(k) + lwrad(k)*(pres(k)/p0)**(R/cp2))
          if ((abs(sd) > R1) .and. (dcond_ls(k) > 0.)) then
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

       ! Either create or evolve ice clouds
       if (create_sgs_ice_clouds(k)) then
          qaiten_create(k) = 0.5/bsi*(bsi+qi_calc)/dt
          qiten_create(k) = qaiten_create(k)*0.5*(bsi+qi_calc) * 0.05 ! AAJ 0.05 to keep ice mass lower
          thiten_create(k) = ((p0/pres(k))**(R/cp2))*lsub*ocp(k)*qiten_create(k)
          ! Create at least 5% cloud fraction
          if ((qa(k) + qaiten_create(k)*dt) < 0.05 .or. (qi(k) + qiten_create(k)*dt) <= R1) then
             qiten_create(k) = 0.
             qaiten_create(k) = 0.
             thiten_create(k) = 0.
          endif
       elseif (evolve_sgs_ice_clouds(k)) then
          ! Contains large-scale forcing and LW cooling (not currently used)
          dcondi_ls(k) = -ali * dqidT(k) * (omega(k)/rho(k)*ocp(k) + lw1d(k)*(pres(k)/p0)**(R/cp2))
          if ((abs(sdi) > R1) .and. (dcondi_ls(k) > 0.)) then
             term1 = ((1.-qa(k))**2*qa(k)**2/qi(k)) + (qa(k)**2*(1.-qa(k))**2/sdi)
             term2 = (1.-qa(k))**2 + qa(k)**2
             gterm = 0.5*term1/term2
             qaiten_evolve(k) = gterm*dcondi_ls(k)
             qiten_evolve(k) = qa(k)*dcondi_ls(k)
             thiten_evolve(k) = ((p0/pres(k))**(R/cp2))*lsub*ocp(k)*qiten_evolve(k)
          else
             qiten_evolve(k) = 0.
             qaiten_evolve(k) = 0.
             thiten_evolve(k) = 0.
          endif
       endif

       ! Erode SGS clouds
       if (erode_sgs_clouds(k)) then
          if ((abs(sd) > R1)) then
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

       if (erode_sgs_ice_clouds(k)) then
          if ((abs(sdi) > R1)) then
             term1 = ((1.-qa(k))**2*qa(k)**2/qi(k)) + (qa(k)**2*(1.-qa(k))**2/sdi)
          else
             term1 = 0.
          endif
          term2 = (1.-qa(k))**2 + qa(k)**2
          gterm = 0.5*term1/term2
          term3 = -3.1*qi_calc/(ali*qvi(k))
          eros_term = -2.25e-5 * exp(term3)
          qiten_erode(k) = min(0., ((qi(k) - qi_calc*qa(k))*eros_term))
          qaiten_erode(k) = min(0., (-gterm*qi_calc*eros_term))
          thiten_erode(k) = ((p0/pres(k))**(R/cp2))*lsub*ocp(k)*qiten_erode(k)
       endif

       ! Sum tendencies
       qaten(k) = qaten_create(k) + qaten_evolve(k) + qaten_erode(k) + &
            qaiten_create(k) + qaiten_evolve(k) + qaiten_erode(k)
       qaten_l(k) = qaten_create(k) + qaten_evolve(k) + qaten_erode(k)
       qaten_i(k) = qaiten_create(k) + qaiten_evolve(k) + qaiten_erode(k)
       qcten(k) = qcten_create(k) + qcten_evolve(k) + qcten_erode(k)
       qiten(k) = qiten_create(k) + qiten_evolve(k) + qiten_erode(k)
       thten(k) = thten_create(k) + thten_evolve(k) + thten_erode(k) + &
            thiten_create(k) + thiten_evolve(k) + thiten_erode(k)
       thten_l(k) = thten_create(k) + thten_evolve(k) + thten_erode(k)
       thten_i(k) = thiten_create(k) + thiten_evolve(k) + thiten_erode(k)

       ! If eroding SGS clouds to less than 5%, or eroding cloud water
       ! to less than R1, then completely remove
       if (((qa(k) + qaten(k)*dt) < 0.05) .and. (qaten(k) < 0.)) then
          qaten(k) = -qa(k)/dt
          qcten(k) = -qc(k)/dt
          qiten(k) = -qi(k)/dt
          thten_l(k) = ((p0/pres(k))**(R/cp2))*lvap(k)*ocp(k)*qcten(k)
          thten_i(k) = ((p0/pres(k))**(R/cp2))*lsub*ocp(k)*qiten(k)
       endif

       if (((qc(k) + qcten(k)*dt) <= R1)) then
!!          qaten(k) = -qa(k)/dt
          qcten(k) = -qc(k)/dt
          thten_l(k) = ((p0/pres(k))**(R/cp2))*lvap(k)*ocp(k)*qcten(k)
       endif

       if (((qi(k) + qiten(k)*dt) <= R1)) then
          qiten(k) = -qi(k)/dt
          thten_i(k) = ((p0/pres(k))**(R/cp2))*lsub*ocp(k)*qiten(k)
          if (((qc(k) + qcten(k)*dt) <= R1)) then
             qaten(k) = -qa(k)/dt
          endif
       endif

       ! Update variables
       qc(k) = qc(k) + qcten(k)*dt
       nc(k) = nc(k) + qcten(k)*dt / (3.14159*1000./6.*((10.e-6)**3.0))
       qi(k) = qi(k) + qiten(k)*dt
       ni(k) = ni(k) + qiten(k)*dt / (3.14159*900./6.*((10.e-6)**3.0))

       qv(k) = qv(k) - qcten(k)*dt - qiten(k)*dt
       qa(k) = qa(k) + qaten(k)*dt
       temp(k) = temp(k) + ((thten_l(k)+thten_i(k)) / ((p0/pres(k))**(R/cp2)))*dt

       if (qc(k) > R1) then
          ! Low-limit check on cloud fraction based on observations
          ! cf = 5.57 * qc [g/kg] ** 0.78
          cf_check = 5.57 * (((qc(k)+qi(k))*1000.)**0.78)
          cf_check = cf_check - (0.5*cf_check)
          cf_check = max(cf_low, min(1., cf_check))

          if (qc(k) >= 2.37e-6) then
             qa(k) = max(qa(k), cf_check)
          endif
       elseif (qi(k) > R1) then
          cf_check = 5.57 * (((10.*qi(k))*1000.)**0.78)
          cf_check = cf_check - (0.5*cf_check)
          cf_check = max(cf_low, min(1., cf_check))

          if (10.*qc(k) >= 2.37e-6) then
             qa(k) = max(qa(k), cf_check)
          endif
       endif
    enddo

  end subroutine tempo_cldfrac_driver

!=================================================================================================================    
end module module_mp_tempo_cldfrac
!=================================================================================================================    
