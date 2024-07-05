module module_mp_thompson_cloud_fraction

  implicit none
  
contains

  subroutine subgrid_scale_ice_clouds(kts, kte, dt, rho, pres, lsub, ocp, qvsi, w1d, &
       no_micro_grid, qa, qv, ri, ni, temp)

    implicit none

    ! Input / Output
    integer, intent(in) :: kts, kte
    real, intent(in) :: rho(:), pres(:), ocp(:), qvsi(:), w1d(:)
    real, intent(in) :: dt, lsub
    logical, intent(in) :: no_micro_grid(:)
    real, intent(inout) :: qa(:), qv(:), ri(:), ni(:), temp(:)

    ! Constants
    real :: con_eps = 0.608
    real :: con_rd = 287.
    real :: grav = 9.8
    real :: cp_air = 1004.
    real :: rvgas = 461.50
    real :: hlv = 2.5e6

    ! Local variables
    integer :: k
    real :: dcond_ls(kts:kte), qa_tend(kts:kte)

    real :: dqsdT(kts:kte), gamm(kts:kte), U(kts:kte), U00(kts:kte), omega(kts:kte)
    real :: dqs_ls(kts:kte), tmp1(kts:kte), tmp2(kts:kte), A_dt(kts:kte), B_dt(kts:kte)
    real :: qa0(kts:kte), qaeq(kts:kte), qa1(kts:kte), qabar(kts:kte), da_ls(kts:kte)
    real :: cooling_terms(kts:kte), eros(kts:kte)
    real :: term1, term2, term3, term4, sink, ratio, pressure_threshold
    logical :: do_subgrid_clouds(kts:kte)

      pressure_threshold = 70000.
      do k = kts, kte
       dcond_ls(k) = 0.
       qa_tend(k) = 0.
       cooling_terms(k) = 0.
       do_subgrid_clouds(k) = .false.

       dqsdT(k) = lsub * qvsi(k) / (rvgas*temp(k)**2)
       gamm(k) = dqsdT(k) * lsub  * ocp(k)
       U(k) = qv(k) / qvsi(k)
       U00(k) = 0.85

       omega(k) = -grav * w1d(k) * rho(k)

       ! Conditions to create sub-grid-scale clouds
       ! Currently hardcoded and arbitary for testing:
       !.. Cloud fraction less than 1
       !.. Temperature between -20 and -40 C
       !.. RH with respect to ice > 85%
       !.. Pressure < 700 hPa
!       if(no_micro_grid(k) .and. temp(k) .le. 253.15 .and. temp(k) .ge. 233.15 &
!            .and. U(k) .ge. U00(k) .and. pres(k) .lt. 70000.) then
!          do_subgrid_clouds(k) = .true.
!       endif

       if(qa(k) .lt. 1.0 .and. temp(k) .le. 253.15 .and. temp(k) .ge. 233.15 &
            .and. U(k) .ge. U00(k) .and. pres(k) .lt. pressure_threshold) then
          do_subgrid_clouds(k) = .true.
       endif

       if (ri(k) .le. 1.e-12) then
          qa(k) = 0.0
       endif
    enddo

    do k = kts, kte
       if (do_subgrid_clouds(k)) then

          if (ri(k) .gt. 1.e-12) then
             qa(k) = max(qa(k), 0.05)
          endif
          ! Currently only have large scale ascent
          ! Need to add convective, radiative, turbulent cooling terms
          cooling_terms(k) = (omega(k)/rho(k)*ocp(k)) * dqsdT(k) * dt

          ! Need net cooling and cloud fraction < 1 for this to work
          if (cooling_terms(k) .lt. 0.0 .and. qa(k) .lt. 1.) then

             term1 = (1.0 + gamm(k)*qa(k))**2.0
             term2 = gamm(k) * (1.0 - qa(k))**2.0
             term3 = 1.0 / qvsi(k) / max(1.0-U(k), 1.e-12)
             term4 = 1.0 + gamm(k)*qa(k)
             tmp1(k) = sqrt(term1 - term2*term3*cooling_terms(k)) - term4
             if (term2*term3 .gt. 1.e-12) then
                tmp1(k) = -1.0 * tmp1(k) / (term2*term3/2.0)
             else
                tmp1(k) = 0.0
             endif
             tmp2(k) = cooling_terms(k) / (1.0 + 0.5*(1.0+qa(k))* gamm(k))
             dqs_ls(k) = min(tmp1(k), tmp2(k))
          else
             dqs_ls(k) = dqs_ls(k)/(1.0 + qa(k)*gamm(k))
          endif

          if (dqs_ls(k) .le. 0.0 .and. qa(k) .le. 1.0) then
             da_ls(k) = -0.5 * (1. - qa(k))*(1. - qa(k)) * dqs_ls(k) * term3
          else
             da_ls(k) = 0.0
          endif

          ! Erosion
          if (ri(k) .gt. 1.e-12) then
             ! 1.e-6 is background erosion (should be scale aware)
             eros(k) = qa(k) * 1.e-6 * dt * qvsi(k) * (1.-U(k)) / ri(k)
          else
             eros(k) = 0.
          endif

          A_dt(k) = 1.0 * da_ls(k) / max(1.0-qa(k), 1.e-12)
          B_dt(k) = eros(k)

          ! Integration with erosion
          if (A_dt(k) .gt. 1.e-12 .or. B_dt(k) .gt. 1.e-12) then
             qa0(k) = qa(k)
             qaeq(k) = A_dt(k) / (A_dt(k) + B_dt(k))
             qa1(k) = qaeq(k) - (qaeq(k) - qa0(k)) * exp(-1.*(A_dt(k)+B_dt(k)))
             qabar(k) = qaeq(k) - ((qa1(k) - qa0(k)) / (A_dt(k) + B_dt(k)))

             ! No erosion
!             qa0(k) = qa(k)
!             qaeq(k) = 1.0
!             qa1(k) = 1.0 - (1.0 - qa(k)) * exp(-1.0*A_dt(k))
!             qabar(k) = 1.0 - (qa1(k) - qa(k)) / A_dt(k)
          else
             qa0(k)   = qa(k)
             qaeq(k)  = qa(k)
             qa1(k)   = qa(k)
             qabar(k) = qa(k)
          endif

          ! Add tendencies for sub-grid clouds
          qa_tend(k) = (qa1(k)-qa0(k)) / dt
          dcond_ls(k) = -1.0 * qabar(k) * dqs_ls(k)
          sink = dcond_ls(k)

          ! Erosions terms (to be implemented later)
          ! dqi/dt = alpha*K*qa*(1-qa)*(qsat-qv)
          ! dqa/dt = -G(-Qc)Qc / (qi-qa*qc)  * dqi/dt
          ! Qc = al * (qt-qsat)
          ! K = 3e-4
          ! alpha = 0.333
          
          ! Limit loss of cloud ice based on available cloud ice
          if (-1.0*dcond_ls(k)*dt .gt. ri(k)) then
             dcond_ls(k) = -ri(k)/dt
          endif

          ! Don't over consume vapor, limit vapor consumption to 90% total vapor
          ratio = 1.0
          if (dcond_ls(k)*dt .gt. 0.9*qv(k)) then
             dcond_ls(k) = 0.9 * qv(k) / dt
             ratio = dcond_ls(k) / sink
          endif
          ! Update the state
          ! Currently assume ice particles nucleate at (arbitrary) 50 microns size 
          qa(k) = qa(k) + qa_tend(k)*dt*ratio
          ri(k) = ri(k) + dcond_ls(k)*dt
          ni(k) = ni(k) + ((dcond_ls(k)*dt) / (3.14159*890./6.*((50.e-6)**3.0)))
          qv(k) = qv(k) - dcond_ls(k)*dt
          temp(k) = temp(k) + lsub*ocp(k)*dcond_ls(k)*dt
          if (ri(k) .lt. 1.e-12) then
             qa(k) = 0.0
             ni(k) = 0.0
          endif
       endif
    enddo
    
  end subroutine subgrid_scale_ice_clouds
  
end module module_mp_thompson_cloud_fraction
