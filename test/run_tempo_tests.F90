program run_tempo_tests
  !! runs tempo tests
  use tests, only : test_tempo_init, test_graupel_sedimentation, &
    test_snow_sedimentation, test_cloud_number_aerosolaware, &
    test_cloud_number_non_aerosolaware, test_cloud_number_ml, &
    test_ml_cloud_effective_radius, test_water_conservation
  use module_mp_tempo_main, only : verify_depletion, depletion_levels_checked, depletion_violations
  
  implicit none

  real, dimension(7) :: sedi_tests = &
    [1., 10., 20., 60., 120., 300., 600.]
  ! maximum |water created or destroyed| / precipitation over the column test;
  ! with sedimentation at the entry density |R/P| is round-off (below 1e-8)
  real, parameter :: water_tolerance = 1.e-6
  integer :: t

  ! tempo init
  call test_tempo_init()
  
  ! ml cloud effective radius
  call test_ml_cloud_effective_radius(dt=20.)
  
  ! test cloud number concentration
  call test_cloud_number_aerosolaware(dt=20.)
  call test_cloud_number_non_aerosolaware(dt=20.)
  call test_cloud_number_ml(dt=20.)

  ! graupel sedimentation
  do t = 1, size(sedi_tests)
    call test_graupel_sedimentation(dt=sedi_tests(t), semi_sedi=.false.)
  enddo
  do t = 1, size(sedi_tests)
    call test_graupel_sedimentation(dt=sedi_tests(t), semi_sedi=.true.)
  enddo

  ! snow sedimentation
  do t = 1, size(sedi_tests)
    call test_snow_sedimentation(dt=sedi_tests(t))
  enddo

  ! during the water conservation runs, also check that no category is over-depleted after
  ! check_over_depletion (net process tendency must not remove more than the category holds)
  verify_depletion = .true.
  depletion_levels_checked = 0
  depletion_violations = 0

  ! column water conservation including precipitation
  do t = 1, size(sedi_tests)
    call test_water_conservation(dt=sedi_tests(t), semi_sedi=.false., tolerance=water_tolerance)
  enddo
  do t = 1, size(sedi_tests)
    call test_water_conservation(dt=sedi_tests(t), semi_sedi=.true., tolerance=water_tolerance)
  enddo

  ! control: sedimentation only (no latent heating, so the density does not change)
  do t = 1, size(sedi_tests)
    call test_water_conservation(dt=sedi_tests(t), semi_sedi=.false., tolerance=water_tolerance, &
      micro_off=.true.)
  enddo
  do t = 1, size(sedi_tests)
    call test_water_conservation(dt=sedi_tests(t), semi_sedi=.true., tolerance=water_tolerance, &
      micro_off=.true.)
  enddo

  ! aerosol-aware: droplet activation, ice nucleation, and aerosol scavenging
  do t = 1, size(sedi_tests)
    call test_water_conservation(dt=sedi_tests(t), semi_sedi=.false., tolerance=water_tolerance, &
      aerosol_aware=.true.)
  enddo
  do t = 1, size(sedi_tests)
    call test_water_conservation(dt=sedi_tests(t), semi_sedi=.true., tolerance=water_tolerance, &
      aerosol_aware=.true.)
  enddo

  ! forced: water vapor added in ascent so cloud water and ice keep forming
  do t = 1, size(sedi_tests)
    call test_water_conservation(dt=sedi_tests(t), semi_sedi=.false., tolerance=water_tolerance, &
      forced=.true.)
  enddo
  do t = 1, size(sedi_tests)
    call test_water_conservation(dt=sedi_tests(t), semi_sedi=.true., tolerance=water_tolerance, &
      forced=.true.)
  enddo
  do t = 1, size(sedi_tests)
    call test_water_conservation(dt=sedi_tests(t), semi_sedi=.false., tolerance=water_tolerance, &
      aerosol_aware=.true., forced=.true.)
  enddo
  do t = 1, size(sedi_tests)
    call test_water_conservation(dt=sedi_tests(t), semi_sedi=.true., tolerance=water_tolerance, &
      aerosol_aware=.true., forced=.true.)
  enddo

  verify_depletion = .false.
  write(*,'(A,I0,A,I0)') 'check_over_depletion --- levels checked: ', depletion_levels_checked, &
    '  over-depleted categories: ', depletion_violations
  if (depletion_levels_checked == 0) error stop 'depletion check did not run'
  if (depletion_violations /= 0) error stop 'a category was over-depleted after check_over_depletion'

end program run_tempo_tests
