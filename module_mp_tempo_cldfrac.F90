! Cloud fraction module for TEMPO Microphysics
!=================================================================================================================
module module_mp_tempo_cldfrac
    use module_mp_tempo_params, only : lvap0, cp2, R, Rv, R1
    use module_mp_tempo_utils, only : rslf
    
#if defined(mpas)
    use mpas_kind_types, only: wp => RKIND, sp => R4KIND, dp => R8KIND
    use mp_radar
#elif defined(standalone)
    use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
    use module_mp_radar
#else
    use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
    use module_mp_radar
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
  subroutine tempo_cldfrac_driver(i,j,kts,kte,dt,temp,pres,rho,w,qa,qv,qc,qi,qaten,qcten,qiten,thten)
!=================================================================================================================
    
    integer, intent(in) :: i, j, kts, kte
    real, intent(in) :: dt
    real, intent(in) :: temp(:), pres(:), rho(:), w(:)
    real, intent(in) :: qv(:), qc(:), qi(:)
    real, intent(in) :: qiten(:)
    real, intent(inout) :: qa(:), qaten(:), qcten(:), thten(:)
    
    integer :: k
    real :: tempc
    real, parameter :: p0 = 100000.
    real, parameter :: critical_rh = 0.85 ! Move params to params module
    real, parameter :: ls_w_limit = 0.1 ! Move params to params module
    real, parameter :: grav = 9.8 ! Move params to params module
    real :: al, bs, sd, qc_calc
    real :: term1, term2, term3, gterm, eros_term

    real, dimension(kts:kte) :: qvs, U, U00, lvap, ocp, dqsdT, omega
    real, dimension(kts:kte) :: cooling_ls
    logical, dimension(kts:kte) :: create_sgs_clouds, erode_sgs_clouds, evolve_sgs_clouds

    do k = kts, kte
       create_sgs_clouds(k) = .false.
       erode_sgs_clouds(k) = .false.
       evolve_sgs_clouds(k) = .false.
       qaten(k) = 0.
       qcten(k) = 0.
       thten(k) = 0.
       cooling_ls(k) = 0.
       qvs(k) = rslf(pres(k), temp(k))
       U(k) = qv(k) / qvs(k)
       U00(k) = critical_rh
       tempc = temp(k) - 273.15       
       lvap(k) = lvap0 + (2106.0 - 4218.0)*tempc
       ocp(k) = 1.0 / (cp2*(1.+0.887*qv(k)))
       omega(k) = -grav * w(k) * rho(k)

       dqsdT(k) = lvap(k) * qvs(k) / (Rv*temp(k)**2)

       ! Limits on cloud fraction based on qc and qi
!       if ((qc(k) <= R1) .and. (qi(k) <= R1)) then
       if ((qc(k) <= R1)) then
          qa(k) = 0.0
       else
          qa(k) = min(qa(k), 1.0)
          qa(k) = max(qa(k), 0.01)
       endif

       ! Make SGS clouds if they do not exist
       ! Only make SGS clouds if cloud water or ice does not already exist
       ! Limit large-scale SGS cloud creation to w < 0.1 m/s
       if((qa(k) < 1.0) .and. (w(k) < ls_w_limit) .and. (U(k) < 1.0) .and. (U(k) > (U00(k) + 0.01)) .and. &
            (qc(k) <= R1)) then
!            (qc(k) <= R1) .and. (qi(k) <= R1)) then
          create_sgs_clouds(k) = .true.
       endif

       ! if(qa(k) < 1.0 .and. qc(k) > R1 .and. (U(k) < 1.) .and. w(k) < ls_w_limit .and. (U(k) > (U00(k) + 0.01))) then
       !    evolve_sgs_clouds(k) = .true.
       ! endif

       ! Erosion condition
!       if(qa(k) < 1.0 .and. ((qc(k) > R1) .or. (qi(k) > R1)) .and. U(k) < 1.) then
       if(qa(k) < 1.0 .and. ((qc(k) > R1)) .and. U(k) < 1.) then
          erode_sgs_clouds(k) = .true.
       endif

       ! Varibles used in cloud fraction scheme
       al = 1. / (1. + dqsdT(k)*lvap(k)*ocp(k))
       bs = al * (1.-U00(k)) * qvs(k)
       sd = al*(qvs(k)-qv(k))
       qc_calc = al*(qv(k)+(qc(k))-qvs(k))

       ! Either create clouds or evolve clouds
       if (create_sgs_clouds(k)) then
          qaten(k) = 0.5/bs*(bs+qc_calc)/dt
          qcten(k) = qaten(k)*0.5*(bs+qc_calc)
          thten(k) = ((p0/pres(k))**(R/cp2))*lvap(k)*ocp(k)*qcten(k)
          ! Create at least 5% cloud fraction
          ! Note that cloud fraction can evolve to be below this value
          if ((qa(k) + qaten(k)*dt) < 0.05) then
             qcten(k) = 0.
             qaten(k) = 0.
             thten(k) = 0.
          endif

!       elseif (evolve_sgs_clouds(k)) then
!          dcond_ls(k) = -al * dqsdT(k) * (omega(k)/rho(k)*ocp(k) + lw1d(k)*(pres(k)/100000.)**0.286 ) * qa(k)
!          if ((abs(sd) > 1.e-12) .and. (dcond_ls(k) > 0.)) then
!             term1 = (1.-qa(k))**2*qa(k)/(rc(k)/rho(k)) + qa(k)**2*(1-qa(k))**2/sd
!             term2 = (1.-qa(k))**2 + qa(k)**2
!             gterm = 0.5*term1/term2
!             qa_tend_ls(k) = gterm/qa(k)*dcond_ls(k)
!          else
!             dcond_ls(k) = 0.
!             qa_tend_ls(k) = 0.
!          endif
       endif

       if (erode_sgs_clouds(k)) then
          if ((abs(sd) > 1.e-12)) then
             term1 = (1.-qa(k))**2*qa(k)**2/(qc(k)) + qa(k)**2*(1-qa(k))**2/sd
          else
             term1 = 0.
          endif

          term2 = (1.-qa(k))**2 + qa(k)**2
          gterm = 0.5*term1/term2

          ! Erosion
          term3 = -3.1*qc_calc/(al*qvs(k))
          eros_term = -2.25e-5 * exp(term3)

          qcten(k) = min(0., (qc(k) - qc_calc*qa(k))*eros_term)
          qaten(k) = min(0., -gterm*qc_calc*eros_term)
          thten(k) = ((p0/pres(k))**(R/cp2))*lvap(k)*ocp(k)*qcten(k)

          if (((qa(k) + qaten(k)*dt) < 0.01) .and. (qaten(k) < 0.)) then
             qaten(k) = -qa(k)/dt
             qcten(k) = -qc(k)/dt
             thten(k) = ((p0/pres(k))**(R/cp2))*lvap(k)*ocp(k)*qcten(k)
          endif

       endif

       !! ALL CF controls need to be here and in the interface
       !! Need to send values from 0.01 to 1 to MP
       !! MP will only set to 1 or 0
    enddo

  end subroutine tempo_cldfrac_driver
  
!=================================================================================================================    
end module module_mp_tempo_cldfrac
!=================================================================================================================    
