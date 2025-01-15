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
  subroutine tempo_cldfrac_driver(i,j,kts,kte,dt,temp,pres,rho,w,qa,qv,qc,qi,nc)
!=================================================================================================================
    
    integer, intent(in) :: i, j, kts, kte
    real, intent(in) :: dt
    real, intent(in) :: pres(:), rho(:), w(:)
    real, intent(inout) :: qa(:), qv(:), qc(:), qi(:), nc(:), temp(:)
    
    integer :: k
    real :: tempc
    real, parameter :: p0 = 100000.
    real, parameter :: ls_w_limit = 0.1
    real, parameter :: grav = 9.8
    real :: al, ali, bs, bsi, sd, sdi, qc_calc
    real :: term1, term2, term3, gterm, eros_term
    real :: cf_check
    real, dimension(kts:kte) :: qaten, qcten, thten, qiten
    real, dimension(kts:kte) :: qvs, qvi, U, Ui, U00, U00i, lvap, ocp, dqsdT, dqidT, omega
    real, dimension(kts:kte) :: qcten_create, qcten_evolve, qcten_erode
    real, dimension(kts:kte) :: qaten_create, qaten_evolve, qaten_erode
    real, dimension(kts:kte) :: thten_create, thten_evolve, thten_erode
    real, dimension(kts:kte) :: cooling_ls, dcond_ls, lw1d
    logical, dimension(kts:kte) :: create_sgs_clouds, erode_sgs_clouds, evolve_sgs_clouds

    do k = kts, kte
       create_sgs_clouds(k) = .false.
       erode_sgs_clouds(k) = .false.
       evolve_sgs_clouds(k) = .false.
       qaten(k) = 0.
       qcten(k) = 0.
       thten(k) = 0.
       qaten_create(k) = 0.
       qcten_create(k) = 0.
       thten_create(k) = 0.
       qaten_evolve(k) = 0.
       qcten_evolve(k) = 0.
       thten_evolve(k) = 0.
       qaten_erode(k) = 0.
       qcten_erode(k) = 0.
       thten_erode(k) = 0.
       dcond_ls(k) = 0.
       lw1d(k) = 0.
       cooling_ls(k) = 0.
       qvs(k) = rslf(pres(k), temp(k))
       qvi(k) = rsif(pres(k), temp(k))
       U(k) = qv(k) / qvs(k)
       Ui(k) = qv(k) / qvi(k)
       U00(k) = critical_rh
       U00i(k) = critical_rh + 0.05
       tempc = temp(k) - 273.15       
       lvap(k) = lvap0 + (2106.0 - 4218.0)*tempc
       ocp(k) = 1.0 / (cp2*(1.+0.887*qv(k)))
       omega(k) = -grav * w(k) * rho(k)
       dqsdT(k) = lvap(k) * qvs(k) / (Rv*temp(k)**2)
       dqidT(k) = lsub * qvi(k) / (Rv*temp(k)**2)

       ! Limit cloud fraction to range of 1 to 100%
       ! if cloud water present
       if ((qc(k) <= R1)) then
          qa(k) = 0.0
       else
          qa(k) = min(qa(k), 1.0)
          qa(k) = max(qa(k), cf_low)
       endif

       ! Create SGS clouds if no cloud water present
       ! Otherwise, evolve SGS clouds from forcing and LW cooling
       if((qa(k) < 1.0) .and. (w(k) < ls_w_limit) .and. (U(k) < 1.0) .and. (U(k) > (U00(k) + 0.01)) .and. tempc > -30.) then
          if (qc(k) <= R1) then
             if (w(k) > 0.) create_sgs_clouds(k) = .true.
          else
             evolve_sgs_clouds(k) = .true.
          endif
       endif

       ! Erode SGS clouds
       if(qa(k) < 1.0 .and. qc(k) > R1 .and. U(k) < 1.) then
          erode_sgs_clouds(k) = .true.
       endif

       ! Varibles used in cloud fraction scheme
       al = 1. / (1. + dqsdT(k)*lvap(k)*ocp(k))
       ali = 1. / (1. + dqidT(k)*lsub*ocp(k))
       bs = al * (1.-U00(k)) * qvs(k)
       bsi = ali * (1.-U00i(k)) * qvi(k)
       sd = al*(qvs(k)-qv(k))
       sdi = ali*(qvi(k)-qv(k))
       qc_calc = al*(qv(k)+qc(k)-qvs(k))

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
          dcond_ls(k) = -al * dqsdT(k) * (omega(k)/rho(k)*ocp(k) + lw1d(k)*(pres(k)/p0)**(R/cp2))
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

       ! Sum tendencies
       qaten(k) = qaten_create(k) + qaten_evolve(k) + qaten_erode(k)
       qcten(k) = qcten_create(k) + qcten_evolve(k) + qcten_erode(k)
       thten(k) = thten_create(k) + thten_evolve(k) + thten_erode(k)

       ! If eroding SGS clouds to less than 5%, or eroding cloud water
       ! to less than R1, then completely remove
       if (((qa(k) + qaten(k)*dt) < 0.05) .and. (qaten(k) < 0.)) then
          qaten(k) = -qa(k)/dt
          qcten(k) = -qc(k)/dt
          thten(k) = ((p0/pres(k))**(R/cp2))*lvap(k)*ocp(k)*qcten(k)
       endif
       if (((qc(k) + qcten(k)*dt) <= R1)) then
          qaten(k) = -qa(k)/dt
          qcten(k) = -qc(k)/dt
          thten(k) = ((p0/pres(k))**(R/cp2))*lvap(k)*ocp(k)*qcten(k)
       endif

       ! Update variables
       qc(k) = qc(k) + qcten(k)*dt
       nc(k) = nc(k) + qcten(k)*dt / (3.14159*1000./6.*((10.e-6)**3.0))
       qv(k) = qv(k) - qcten(k)*dt
       qa(k) = qa(k) + qaten(k)*dt
       temp(k) = temp(k) + (thten(k) / ((p0/pres(k))**(R/cp2)))*dt

       ! Low-limit check on cloud fraction based on observations
       ! cf = 5.57 * qc [g/kg] ** 0.78
       cf_check = 5.57 * ((qc(k)*1000.)**0.78)
       cf_check = cf_check - (0.5*cf_check)
       cf_check = max(cf_low, min(1., cf_check))

       if (qc(k) >= 2.37e-6) then
          qa(k) = max(qa(k), cf_check)
       endif
    enddo

  end subroutine tempo_cldfrac_driver
  
!=================================================================================================================    
end module module_mp_tempo_cldfrac
!=================================================================================================================    
