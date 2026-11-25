module tests
  !! TEMPO tests
  use module_mp_tempo_cfgs, only : ty_tempo_cfgs
  use module_mp_tempo_driver, only : tempo_init, tempo_run, ty_tempo_driver_diags
  implicit none
  private

  public :: test_tempo_init, test_graupel_sedimentation, test_snow_sedimentation, &
    test_cloud_number_aerosolaware, test_cloud_number_non_aerosolaware, &
    test_cloud_number_ml, test_ml_cloud_effective_radius, test_water_conservation, &
    test_cloud_fraction_conservation

  type(ty_tempo_cfgs) :: tempo_cfgs
  contains

  subroutine test_tempo_init()
    !! test tempo initialization procedure
    !! use ml_for_bl_nc_flag = .true. to initialize ml data
    !! which is used for a few tests
    !! test specific flags are then set in each test
    call tempo_init(ml_for_bl_nc_flag = .true., tempo_cfgs=tempo_cfgs)
  end subroutine test_tempo_init


 subroutine test_graupel_sedimentation(dt, semi_sedi)
    !! test graupel sedimentation
    use module_mp_tempo_params, only : wp, sp, dp

    integer :: itimestep
    real(wp), intent(in) :: dt
    logical, intent(in) :: semi_sedi
    integer, parameter :: nz = 59
    integer, parameter :: ids = 1, ide = 1, ims = 1, ime = 1, its = 1, ite = 1
    integer, parameter :: jds = 1, jde = 1, jms = 1, jme = 1, jts = 1, jte = 1
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(its:ite, kts:kte, jts:jte) :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr, ng, qb
    real(wp), dimension(nz) :: klevs_in, qv_in, qc_in, qr_in, qi_in, qs_in, &
      qg_in, ni_in, nr_in, nc_in, nwfa_in, nifa_in, theta_in, ng_in, volg_in, &
      pressure_in, w_in, dz_in
    real(wp) :: precip_sum
    character(len=20) :: dt_string, tt_string
    character(len=20) :: semi_sedi_string
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    type(ty_tempo_cfgs) :: tempo_cfgs
    integer :: io1, io2, io3, io4, k, total_timesteps
    integer, parameter :: integration_time = 1200

    write(dt_string, '(I7)') int(dt)
    
    if (semi_sedi) then
      write(semi_sedi_string, '(A)') '_semi_sedi'
    else
      write(semi_sedi_string, '(A)') ''
    endif

    ! read input file
    io1 = 11
    open (io1, file='../test/data/mpas_59lev_test.txt', status='old')
    read(io1,*) ! header
    do k = 1, nz
      read(io1,*) klevs_in(k), qv_in(k), qc_in(k), qr_in(k), qi_in(k), qs_in(k), &
          qg_in(k), ni_in(k), nr_in(k), nc_in(k), nwfa_in(k), nifa_in(k), &
          theta_in(k), ng_in(k), volg_in(k), pressure_in(k), w_in(k), dz_in(k)
    end do
    close(io1)

    qv(1,:,1) = qv_in
    t(1,:,1) = theta_in * (pressure_in/100000.)**0.286
    p(1,:,1) = pressure_in
    dz(1,:,1) = dz_in
    w = 0._wp
    qc = 0._wp
    qr = 0._wp
    nr = 0._wp
    qi = 0._wp
    qs = 0._wp
    ni = 0._wp
    qg(1,:,1) = qg_in
    ng(1,:,1) = ng_in
    qb(1,:,1) = volg_in * 1000. ! convert meters^3 -> liters

    ! set configs
    tempo_cfgs%turn_off_micro_flag = .true.
    tempo_cfgs%graupel_med_vol_diam_flag = .true.
    tempo_cfgs%semi_sedi_flag = semi_sedi

    total_timesteps = int(integration_time/dt)
    precip_sum = 0._wp
    open(newunit=io2, file="graupel_precip_dt_"//trim(adjustl(dt_string))//""//trim(semi_sedi_string)//".txt", &
      status="new", action="write")

    do itimestep = 1, total_timesteps
      call  tempo_run(tempo_cfgs=tempo_cfgs, itimestep=itimestep, dt=dt, &
                      ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                      jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                      kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                      qc=qc, qr=qr, qi=qi, qs=qs, qg=qg,  ni=ni, nr=nr, ng=ng, qb=qb, &
                      tempo_diags=tempo_driver_diags)
      precip_sum = precip_sum + tempo_driver_diags%graupel_liquid_equiv_precip(1,1)
      write(io2,'(I7, 1E12.4)') int(itimestep*dt), precip_sum

      if (dt == 1. .and. .not. semi_sedi .and. itimestep == 1) then
        open(newunit=io3, file="graupel_sedi_init.txt", status="new", action="write")
        write(io3, '(5A)') 'k ', 'mass ', 'number ', 'density ' , 'mvd'
        do k = 1, nz
          if (qb(1,k,1) > 0._wp) then
            write(io3,'(I5, 4E12.4)') k, qg(1,k,1), ng(1,k,1), 1000.*qg(1,k,1)/qb(1,k,1), tempo_driver_diags%graupel_med_vol_diam(1,k,1)
          else
            write(io3,'(I5, 4E12.4)') k, qg(1,k,1), ng(1,k,1), 0._wp, tempo_driver_diags%graupel_med_vol_diam(1,k,1)
          endif 
        enddo 
      endif 
    enddo

    write(tt_string, '(I7)') integration_time
    open(newunit=io4, file="graupel_sedi_dt_"//trim(adjustl(dt_string))//""//trim(semi_sedi_string)//&
      "_runtime"//trim(adjustl(tt_string))//".txt", status="new", action="write")
      write(io4, '(5A)') 'k ', 'mass ', 'number ', 'density ' , 'mvd'
    do k = 1, nz
      if (qb(1,k,1) > 0._wp) then
        write(io4,'(I5, 4E12.4)') k, qg(1,k,1), ng(1,k,1), 1000.*qg(1,k,1)/qb(1,k,1), &
          tempo_driver_diags%graupel_med_vol_diam(1,k,1)
      else
        write(io4,'(I5, 4E12.4)') k, qg(1,k,1), ng(1,k,1), 0._wp, &
          tempo_driver_diags%graupel_med_vol_diam(1,k,1)
      endif 
    enddo 
  end subroutine test_graupel_sedimentation


  subroutine test_snow_sedimentation(dt)
    use module_mp_tempo_params, only : wp, sp, dp

    integer :: itimestep
    real(wp), intent(in) :: dt
    integer, parameter :: nz = 59
    integer, parameter :: ids = 1, ide = 1, ims = 1, ime = 1, its = 1, ite = 1
    integer, parameter :: jds = 1, jde = 1, jms = 1, jme = 1, jts = 1, jte = 1
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(its:ite, kts:kte, jts:jte) :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr, ng, qb
    real(wp), dimension(nz) :: klevs_in, qv_in, qc_in, qr_in, qi_in, qs_in, &
      qg_in, ni_in, nr_in, nc_in, nwfa_in, nifa_in, theta_in, ng_in, volg_in, &
      pressure_in, w_in, dz_in
    real(wp) :: precip_sum
    character(len=20) :: dt_string, tt_string
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    type(ty_tempo_cfgs) :: tempo_cfgs
    integer :: io1, io2, io3, io4, k, total_timesteps
    integer, parameter :: integration_time = 1200

    write(dt_string, '(I7)') int(dt)

    ! read input file
    io1 = 11
    open (io1, file='../test/data/mpas_59lev_test.txt', status='old')
    read(io1,*) ! header
    do k = 1, nz
      read(io1,*) klevs_in(k), qv_in(k), qc_in(k), qr_in(k), qi_in(k), qs_in(k), &
          qg_in(k), ni_in(k), nr_in(k), nc_in(k), nwfa_in(k), nifa_in(k), &
          theta_in(k), ng_in(k), volg_in(k), pressure_in(k), w_in(k), dz_in(k)
    end do
    close(io1)

    qv(1,:,1) = qv_in
    t(1,:,1) = theta_in * (pressure_in/100000.)**0.286
    p(1,:,1) = pressure_in
    dz(1,:,1) = dz_in
    w = 0._wp
    qc = 0._wp
    qr = 0._wp
    nr = 0._wp
    qi = 0._wp
    qs(1,:,1) = qs_in
    ni = 0._wp
    qg = 0._wp
    ng = 0._wp
    qb = 0._wp

    ! set configs
    tempo_cfgs%turn_off_micro_flag= .true.

    total_timesteps = int(integration_time/dt)
    precip_sum = 0._wp
    open(newunit=io2, file="snow_precip_dt_"//trim(adjustl(dt_string))//".txt", &
      status="new", action="write")

    do itimestep = 1, total_timesteps
      call  tempo_run(tempo_cfgs=tempo_cfgs, itimestep=itimestep, dt=dt, &
                      ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                      jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                      kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                      qc=qc, qr=qr, qi=qi, qs=qs, qg=qg,  ni=ni, nr=nr, ng=ng, qb=qb, &
                      tempo_diags=tempo_driver_diags)
      precip_sum = precip_sum + tempo_driver_diags%snow_liquid_equiv_precip(1,1)
      write(io2,'(I7, 1E12.4)') int(itimestep*dt), precip_sum

      if (dt == 1. .and. itimestep == 1) then
        open(newunit=io3, file="snow_sedi_init.txt", status="new", action="write")
        do k = 1, nz
          write(io3,'(I5, 4E12.4)') k, qs(1,k,1)
        enddo 
      endif 
    enddo

    write(tt_string, '(I7)') integration_time
    open(newunit=io4, file="snow_sedi_dt_"//trim(adjustl(dt_string))// &
      "_runtime"//trim(adjustl(tt_string))//".txt", status="new", action="write")
    do k = 1, nz
      write(io4,'(I5, 4E12.4)') k, qs(1,k,1)
    enddo 
  end subroutine test_snow_sedimentation


  subroutine test_cloud_number_aerosolaware(dt)
    use module_mp_tempo_params, only : wp, sp, dp

    integer :: itimestep
    real(wp), intent(in) :: dt
    integer, parameter :: nz = 59
    integer, parameter :: ids = 1, ide = 1, ims = 1, ime = 1, its = 1, ite = 1
    integer, parameter :: jds = 1, jde = 1, jms = 1, jme = 1, jts = 1, jte = 1
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(its:ite, kts:kte, jts:jte) :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr, ng, qb, nc
    real(wp), dimension(nz) :: klevs_in, qv_in, qc_in, qr_in, qi_in, qs_in, &
      qg_in, ni_in, nr_in, nc_in, nwfa_in, nifa_in, theta_in, ng_in, volg_in, &
      pressure_in, w_in, dz_in
    real(wp) :: precip_sum
    character(len=20) :: dt_string, tt_string
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    type(ty_tempo_cfgs) :: tempo_cfgs
    integer :: io1, io2, io3, io4, k, total_timesteps
    integer, parameter :: integration_time = 1200

    write(dt_string, '(I7)') int(dt)

    ! read input file
    io1 = 11
    open (io1, file='../test/data/mpas_59lev_test.txt', status='old')
    read(io1,*) ! header
    do k = 1, nz
      read(io1,*) klevs_in(k), qv_in(k), qc_in(k), qr_in(k), qi_in(k), qs_in(k), &
          qg_in(k), ni_in(k), nr_in(k), nc_in(k), nwfa_in(k), nifa_in(k), &
          theta_in(k), ng_in(k), volg_in(k), pressure_in(k), w_in(k), dz_in(k)
    end do
    close(io1)

    qv(1,:,1) = qv_in
    t(1,:,1) = theta_in * (pressure_in/100000.)**0.286
    p(1,:,1) = pressure_in
    dz(1,:,1) = dz_in
    w(1,:,1) = w_in
    qc(1,:,1) = qc_in
    nc(1,:,1) = nc_in
    qr(1,:,1) = qr_in
    nr(1,:,1) = nr_in
    qi(1,:,1) =  qi_in
    qs(1,:,1) = qs_in
    ni(1,:,1) = ni_in
    qg = 0._wp
    ng = 0._wp
    qb = 0._wp

    ! set configs
    tempo_cfgs%turn_off_micro_flag= .true.

    total_timesteps = int(integration_time/dt)

    do itimestep = 1, total_timesteps
      call  tempo_run(tempo_cfgs=tempo_cfgs, itimestep=itimestep, dt=dt, &
                      ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                      jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                      kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                      qc=qc, nc=nc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
                      tempo_diags=tempo_driver_diags)

      if (itimestep == 1) then
        open(newunit=io3, file="cloud_number_init.txt", status="new", action="write")
        do k = 1, nz
          write(io3,'(I5, 4E12.4)') k, qc(1,k,1), nc(1,k,1)
        enddo 
      endif 
    enddo

    write(tt_string, '(I7)') integration_time
    open(newunit=io4, file="cloud_number_dt_"//trim(adjustl(dt_string))// &
      "_runtime"//trim(adjustl(tt_string))//".txt", status="new", action="write")
    do k = 1, nz
      write(io4,'(I5, 4E12.4)') k, qc(1,k,1), nc(1,k,1)
    enddo 
  end subroutine test_cloud_number_aerosolaware


  subroutine test_cloud_number_non_aerosolaware(dt)
    use module_mp_tempo_params, only : wp, sp, dp

    integer :: itimestep
    real(wp), intent(in) :: dt
    integer, parameter :: nz = 59
    integer, parameter :: ids = 1, ide = 1, ims = 1, ime = 1, its = 1, ite = 1
    integer, parameter :: jds = 1, jde = 1, jms = 1, jme = 1, jts = 1, jte = 1
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(its:ite, kts:kte, jts:jte) :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr, ng, qb, nc
    real(wp), dimension(nz) :: klevs_in, qv_in, qc_in, qr_in, qi_in, qs_in, &
      qg_in, ni_in, nr_in, nc_in, nwfa_in, nifa_in, theta_in, ng_in, volg_in, &
      pressure_in, w_in, dz_in
    real(wp) :: precip_sum
    character(len=20) :: dt_string, tt_string
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    type(ty_tempo_cfgs) :: tempo_cfgs
    integer :: io1, io2, io3, io4, k, total_timesteps
    integer, parameter :: integration_time = 1200

    write(dt_string, '(I7)') int(dt)

    ! read input file
    io1 = 11
    open (io1, file='../test/data/mpas_59lev_test.txt', status='old')
    read(io1,*) ! header
    do k = 1, nz
      read(io1,*) klevs_in(k), qv_in(k), qc_in(k), qr_in(k), qi_in(k), qs_in(k), &
          qg_in(k), ni_in(k), nr_in(k), nc_in(k), nwfa_in(k), nifa_in(k), &
          theta_in(k), ng_in(k), volg_in(k), pressure_in(k), w_in(k), dz_in(k)
    end do
    close(io1)

    qv(1,:,1) = qv_in
    t(1,:,1) = theta_in * (pressure_in/100000.)**0.286
    p(1,:,1) = pressure_in
    dz(1,:,1) = dz_in
    w(1,:,1) = w_in
    qc(1,:,1) = qc_in
    nc(1,:,1) = 0._wp
    qr(1,:,1) = qr_in
    nr(1,:,1) = nr_in
    qi(1,:,1) =  qi_in
    qs(1,:,1) = qs_in
    ni(1,:,1) = ni_in
    qg = 0._wp
    ng = 0._wp
    qb = 0._wp

    ! set configs
    tempo_cfgs%turn_off_micro_flag= .true.
    tempo_cfgs%cloud_number_mixing_ratio_flag= .true.
    tempo_cfgs%ml_for_nc_flag = .false.

    total_timesteps = int(integration_time/dt)

    do itimestep = 1, total_timesteps
      call  tempo_run(tempo_cfgs=tempo_cfgs, itimestep=itimestep, dt=dt, &
                      ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                      jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                      kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                      qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
                      tempo_diags=tempo_driver_diags)

      if (itimestep == 1) then
        open(newunit=io3, file="cloud_number_constant_init.txt", status="new", action="write")
        do k = 1, nz
          write(io3,'(I5, 4E12.4)') k, qc(1,k,1), tempo_driver_diags%cloud_number_mixing_ratio(1,k,1)
        enddo 
      endif 
    enddo

    write(tt_string, '(I7)') integration_time
    open(newunit=io4, file="cloud_number_constant_dt_"//trim(adjustl(dt_string))// &
      "_runtime"//trim(adjustl(tt_string))//".txt", status="new", action="write")
    do k = 1, nz
      write(io4,'(I5, 4E12.4)') k, qc(1,k,1), tempo_driver_diags%cloud_number_mixing_ratio(1,k,1)
    enddo 
  end subroutine test_cloud_number_non_aerosolaware


  subroutine test_cloud_number_ml(dt)
    use module_mp_tempo_params, only : wp, sp, dp

    integer :: itimestep
    real(wp), intent(in) :: dt
    integer, parameter :: nz = 59
    integer, parameter :: ids = 1, ide = 1, ims = 1, ime = 1, its = 1, ite = 1
    integer, parameter :: jds = 1, jde = 1, jms = 1, jme = 1, jts = 1, jte = 1
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(its:ite, kts:kte, jts:jte) :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr, ng, qb, nc
    real(wp), dimension(nz) :: klevs_in, qv_in, qc_in, qr_in, qi_in, qs_in, &
      qg_in, ni_in, nr_in, nc_in, nwfa_in, nifa_in, theta_in, ng_in, volg_in, &
      pressure_in, w_in, dz_in
    real(wp) :: precip_sum
    character(len=20) :: dt_string, tt_string
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    type(ty_tempo_cfgs) :: tempo_cfgs
    integer :: io1, io2, io3, io4, k, total_timesteps
    integer, parameter :: integration_time = 1200

    write(dt_string, '(I7)') int(dt)

    ! read input file
    io1 = 11
    open (io1, file='../test/data/mpas_59lev_test.txt', status='old')
    read(io1,*) ! header
    do k = 1, nz
      read(io1,*) klevs_in(k), qv_in(k), qc_in(k), qr_in(k), qi_in(k), qs_in(k), &
          qg_in(k), ni_in(k), nr_in(k), nc_in(k), nwfa_in(k), nifa_in(k), &
          theta_in(k), ng_in(k), volg_in(k), pressure_in(k), w_in(k), dz_in(k)
    end do
    close(io1)

    qv(1,:,1) = qv_in
    t(1,:,1) = theta_in * (pressure_in/100000.)**0.286
    p(1,:,1) = pressure_in
    dz(1,:,1) = dz_in
    w(1,:,1) = w_in
    qc(1,:,1) = qc_in
    nc(1,:,1) = 0._wp
    qr(1,:,1) = qr_in
    nr(1,:,1) = nr_in
    qi(1,:,1) =  qi_in
    qs(1,:,1) = qs_in
    ni(1,:,1) = ni_in
    qg = 0._wp
    ng = 0._wp
    qb = 0._wp

    ! set configs
    tempo_cfgs%turn_off_micro_flag = .true.
    tempo_cfgs%cloud_number_mixing_ratio_flag = .true.
    tempo_cfgs%ml_for_nc_flag = .true.

    total_timesteps = int(integration_time/dt)

    do itimestep = 1, total_timesteps
      call  tempo_run(tempo_cfgs=tempo_cfgs, itimestep=itimestep, dt=dt, &
                      ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                      jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                      kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                      qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
                      tempo_diags=tempo_driver_diags)

      if (itimestep == 1) then
        open(newunit=io3, file="cloud_number_ml_init.txt", status="new", action="write")
        do k = 1, nz
          write(io3,'(I5, 4E12.4)') k, qc(1,k,1), tempo_driver_diags%cloud_number_mixing_ratio(1,k,1)
        enddo 
      endif 
    enddo

    write(tt_string, '(I7)') integration_time
    open(newunit=io4, file="cloud_number_ml_dt_"//trim(adjustl(dt_string))// &
      "_runtime"//trim(adjustl(tt_string))//".txt", status="new", action="write")
    do k = 1, nz
      write(io4,'(I5, 4E12.4)') k, qc(1,k,1), tempo_driver_diags%cloud_number_mixing_ratio(1,k,1)
    enddo 
  end subroutine test_cloud_number_ml


  subroutine test_ml_cloud_effective_radius(dt)
    use module_mp_tempo_params, only : wp, sp, dp

    integer :: itimestep
    real(wp), intent(in) :: dt
    integer, parameter :: nz = 59
    integer, parameter :: ids = 1, ide = 1, ims = 1, ime = 1, its = 1, ite = 1
    integer, parameter :: jds = 1, jde = 1, jms = 1, jme = 1, jts = 1, jte = 1
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(its:ite, kts:kte, jts:jte) :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr, ng, qb, nc, qc_bl, qcfrac_bl
    real(wp), dimension(nz) :: klevs_in, qv_in, qc_in, qr_in, qi_in, qs_in, &
      qg_in, ni_in, nr_in, nc_in, nwfa_in, nifa_in, theta_in, ng_in, volg_in, &
      pressure_in, w_in, dz_in
    real(wp) :: precip_sum
    character(len=20) :: dt_string, tt_string
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    type(ty_tempo_cfgs) :: tempo_cfgs
    integer :: io1, io2, io3, io4, k, total_timesteps
    integer, parameter :: integration_time = 1200

    write(dt_string, '(I7)') int(dt)

    ! read input file
    io1 = 11
    open (io1, file='../test/data/mpas_59lev_test.txt', status='old')
    read(io1,*) ! header
    do k = 1, nz
      read(io1,*) klevs_in(k), qv_in(k), qc_in(k), qr_in(k), qi_in(k), qs_in(k), &
          qg_in(k), ni_in(k), nr_in(k), nc_in(k), nwfa_in(k), nifa_in(k), &
          theta_in(k), ng_in(k), volg_in(k), pressure_in(k), w_in(k), dz_in(k)
    end do
    close(io1)

    qv(1,:,1) = qv_in
    t(1,:,1) = theta_in * (pressure_in/100000.)**0.286
    p(1,:,1) = pressure_in
    dz(1,:,1) = dz_in
    w(1,:,1) = w_in
    qc(1,:,1) = qc_in
    nc(1,:,1) = nc_in
    qr(1,:,1) = qr_in
    nr(1,:,1) = nr_in
    qi(1,:,1) =  qi_in
    qs(1,:,1) = qs_in
    ni(1,:,1) = ni_in

    qc_bl(1,:,1) = 0._wp
    qcfrac_bl(1,:,1) = 0._wp
    do k = 1, nz
      if (k < 20) then
        qc_bl(1,k,1) = qc_in(k+15)*0.2_wp
      endif 
      if(qc_bl(1,k,1) > 0._wp) then
        qcfrac_bl(1,k,1) = 0.35_wp
      endif 
    enddo 
    qg = 0._wp
    ng = 0._wp
    qb = 0._wp

    ! set configs
    tempo_cfgs%turn_off_micro_flag = .true.
    tempo_cfgs%ml_for_bl_nc_flag = .true.
    tempo_cfgs%cloud_number_mixing_ratio_flag = .true.

    total_timesteps = int(integration_time/dt)

    do itimestep = 1, total_timesteps
      call  tempo_run(tempo_cfgs=tempo_cfgs, itimestep=itimestep, dt=dt, &
                      ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                      jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                      kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                      qc_bl=qc_bl, qcfrac_bl=qcfrac_bl, &
                      qc=qc, nc=nc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
                      tempo_diags=tempo_driver_diags)

      if (itimestep == 1) then
        open(newunit=io3, file="cloud_re_init.txt", status="new", action="write")
        do k = 1, nz
          write(io3,'(I5, 5E12.4)') k, qc(1,k,1), qc_bl(1,k,1), qcfrac_bl(1,k,1), &
            tempo_driver_diags%cloud_number_mixing_ratio(1,k,1), tempo_driver_diags%re_cloud(1,k,1)*1.e6_wp
        enddo 
      endif 
    enddo

    write(tt_string, '(I7)') integration_time
    open(newunit=io4, file="cloud_re_dt_"//trim(adjustl(dt_string))// &
      "_runtime"//trim(adjustl(tt_string))//".txt", status="new", action="write")
    do k = 1, nz
      write(io4,'(I5, 5E12.4)') k, qc(1,k,1), qc_bl(1,k,1), qcfrac_bl(1,k,1), &
        tempo_driver_diags%cloud_number_mixing_ratio(1,k,1), tempo_driver_diags%re_cloud(1,k,1)*1.e6_wp
    enddo 
  end subroutine test_ml_cloud_effective_radius


  subroutine test_water_conservation(dt, semi_sedi, tolerance, micro_off, aerosol_aware, forced, &
    residual_ratio)
    !! column water budget with microphysics and precipitation
    !!
    !! each step, total water W = sum(m * (qv+qc+qr+qi+qs+qg)) is computed before and
    !! after tempo_run using the same layer mass m = rho*dz, where rho is the density
    !! TEMPO computes from the state on entry to the step. With P the surface precipitation
    !! (rain, cloud, ice, snow, and graupel/hail) in kg m^-2 during the step, the water the
    !! scheme creates (+) or destroys (-) is R = (W_after - W_before) + P.
    !!
    !! the test reports the accumulated R as a fraction of the accumulated P and stops with
    !! an error if |R/P| exceeds tolerance
    !!
    !! options:
    !!   micro_off = .true. turns off all microphysical processes so only sedimentation
    !!     acts; with no latent heating the density does not change, so R is round-off
    !!     and a failure points at sedimentation itself
    !!   aerosol_aware = .true. passes water- and ice-friendly aerosols (nwfa, nifa) so
    !!     droplet activation, ice nucleation, and aerosol scavenging are used
    !!   forced = .true. adds water vapor where w > 0 before each step, a stand-in for
    !!     the moistening by ascent that a host model's dynamics would provide, at
    !!     dqv/dt = forcing_gamma * min(w, forcing_w_max) throughout the run, so cloud
    !!     water and ice keep forming. The forcing is applied outside tempo_run,
    !!     before W is computed, so the budget for the step is unchanged; the injected
    !!     water F (kg m^-2) is reported
    use module_mp_tempo_params, only : wp, sp, dp, roverrv, rdry, min_qv

    real(wp), intent(in) :: dt
    logical, intent(in) :: semi_sedi
    real(wp), intent(in) :: tolerance
    logical, intent(in), optional :: micro_off, aerosol_aware, forced
    real(wp), intent(out), optional :: residual_ratio
    integer :: itimestep
    integer, parameter :: nz = 59
    integer, parameter :: ids = 1, ide = 1, ims = 1, ime = 1, its = 1, ite = 1
    integer, parameter :: jds = 1, jde = 1, jms = 1, jme = 1, jts = 1, jte = 1
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(its:ite, kts:kte, jts:jte) :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr, ng, qb, nc, nwfa, nifa
    real(wp), dimension(nz) :: klevs_in, qv_in, qc_in, qr_in, qi_in, qs_in, &
      qg_in, ni_in, nr_in, nc_in, nwfa_in, nifa_in, theta_in, ng_in, volg_in, &
      pressure_in, w_in, dz_in
    real(wp), dimension(nz) :: layer_mass, forcing_step
    real(dp) :: water_before, water_after, precip_step, residual_step, &
      precip_sum, residual_sum, max_rel_step, ratio, forcing_sum, cloud_end
    real(wp), parameter :: forcing_gamma = 1.e-6_wp !! vapor added per meter of ascent \([kg\; kg^{-1}\; m^{-1}]\)
    real(wp), parameter :: forcing_w_max = 4._wp !! vertical velocity cap for the forcing \([m\; s^{-1}]\)
    character(len=20) :: dt_string
    character(len=20) :: semi_sedi_string
    character(len=40) :: case_string
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    type(ty_tempo_cfgs) :: tempo_cfgs
    integer :: io1, io2, k, total_timesteps
    integer, parameter :: integration_time = 1200
    logical :: no_micro, aerosols, add_forcing

    add_forcing = .false.
    if (present(forced)) add_forcing = forced
    no_micro = .false.
    if (present(micro_off)) no_micro = micro_off
    aerosols = .false.
    if (present(aerosol_aware)) aerosols = aerosol_aware

    write(dt_string, '(I7)') int(dt)
    if (semi_sedi) then
      write(semi_sedi_string, '(A)') '_semi_sedi'
    else
      write(semi_sedi_string, '(A)') ''
    endif
    case_string = ''
    if (no_micro) case_string = trim(case_string)//'_micro_off'
    if (aerosols) case_string = trim(case_string)//'_aerosol_aware'
    if (add_forcing) case_string = trim(case_string)//'_forced'

    ! read input file
    open(newunit=io1, file='../test/data/mpas_59lev_test.txt', status='old')
    read(io1,*) ! header
    do k = 1, nz
      read(io1,*) klevs_in(k), qv_in(k), qc_in(k), qr_in(k), qi_in(k), qs_in(k), &
          qg_in(k), ni_in(k), nr_in(k), nc_in(k), nwfa_in(k), nifa_in(k), &
          theta_in(k), ng_in(k), volg_in(k), pressure_in(k), w_in(k), dz_in(k)
    end do
    close(io1)

    ! all hydrometeors present and all processes on
    qv(1,:,1) = qv_in
    t(1,:,1) = theta_in * (pressure_in/100000.)**0.286
    p(1,:,1) = pressure_in
    dz(1,:,1) = dz_in
    w(1,:,1) = w_in
    qc(1,:,1) = qc_in
    nc(1,:,1) = nc_in
    qr(1,:,1) = qr_in
    nr(1,:,1) = nr_in
    qi(1,:,1) = qi_in
    ni(1,:,1) = ni_in
    qs(1,:,1) = qs_in
    qg(1,:,1) = qg_in
    ng(1,:,1) = ng_in
    qb(1,:,1) = volg_in * 1000. ! convert meters^3 -> liters
    nwfa(1,:,1) = nwfa_in
    nifa(1,:,1) = nifa_in

    ! set configs
    tempo_cfgs%turn_off_micro_flag = no_micro
    tempo_cfgs%semi_sedi_flag = semi_sedi

    total_timesteps = int(integration_time/dt)
    precip_sum = 0._dp
    residual_sum = 0._dp
    max_rel_step = 0._dp
    forcing_sum = 0._dp
    forcing_step = 0._wp
    open(newunit=io2, file="water_budget_dt_"//trim(adjustl(dt_string))//trim(semi_sedi_string)// &
      trim(case_string)//".txt", &
      status="new", action="write")
    write(io2, '(A)') 'time water_column precip_sum residual_sum forcing_sum'

    do itimestep = 1, total_timesteps
      ! forcing: add water vapor in ascent, outside tempo_run
      if (add_forcing) then
        do k = 1, nz
          forcing_step(k) = forcing_gamma*min(max(w(1,k,1), 0._wp), forcing_w_max)*dt
          qv(1,k,1) = qv(1,k,1) + forcing_step(k)
        enddo
      endif

      ! layer mass from the density TEMPO computes on entry to the step
      do k = 1, nz
        layer_mass(k) = roverrv*p(1,k,1) / (rdry*t(1,k,1)*(max(min_qv, qv(1,k,1))+roverrv)) * dz(1,k,1)
      enddo
      water_before = column_water(layer_mass, qv(1,:,1), qc(1,:,1), qr(1,:,1), &
        qi(1,:,1), qs(1,:,1), qg(1,:,1))
      forcing_sum = forcing_sum + dot_product(real(layer_mass, kind=dp), real(forcing_step, kind=dp))

      if (aerosols) then
        call  tempo_run(tempo_cfgs=tempo_cfgs, itimestep=itimestep, dt=dt, &
                        ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                        jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                        kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                        t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                        qc=qc, nc=nc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, ng=ng, qb=qb, &
                        nwfa=nwfa, nifa=nifa, tempo_diags=tempo_driver_diags)
      else
        call  tempo_run(tempo_cfgs=tempo_cfgs, itimestep=itimestep, dt=dt, &
                        ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                        jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                        kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                        t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                        qc=qc, nc=nc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, ng=ng, qb=qb, &
                        tempo_diags=tempo_driver_diags)
      endif

      water_after = column_water(layer_mass, qv(1,:,1), qc(1,:,1), qr(1,:,1), &
        qi(1,:,1), qs(1,:,1), qg(1,:,1))

      ! graupel_liquid_equiv_precip includes hail
      precip_step = tempo_driver_diags%rain_precip(1,1) + tempo_driver_diags%cloud_precip(1,1) + &
        tempo_driver_diags%ice_liquid_equiv_precip(1,1) + tempo_driver_diags%snow_liquid_equiv_precip(1,1) + &
        tempo_driver_diags%graupel_liquid_equiv_precip(1,1)

      residual_step = (water_after - water_before) + precip_step
      precip_sum = precip_sum + precip_step
      residual_sum = residual_sum + residual_step
      max_rel_step = max(max_rel_step, abs(residual_step)/water_before)
      write(io2, '(I7, 4E16.8)') int(itimestep*dt), water_after, precip_sum, residual_sum, forcing_sum
    enddo
    close(io2)

    ! cloud water + ice left at the end
    cloud_end = dot_product(real(layer_mass, kind=dp), real(qc(1,:,1) + qi(1,:,1), kind=dp))

    ratio = 0._dp
    if (precip_sum > 0._dp) ratio = residual_sum / precip_sum
    if (present(residual_ratio)) residual_ratio = real(ratio, kind=wp)

    write(*,'(A,I5,A,L1,A,L1,A,L1,A,L1,A,ES11.4,A,ES11.4,A,ES11.4,A,ES11.4,A,ES10.2,A,ES10.2)') &
      'test_water_conservation() --- dt = ', int(dt), ' semi_sedi = ', semi_sedi, &
      ' micro_off = ', no_micro, ' aerosol_aware = ', aerosols, ' forced = ', add_forcing, &
      ' F = ', forcing_sum, ' P = ', precip_sum, ' qc+qi end = ', cloud_end, &
      ' kg m-2  R = ', residual_sum, ' kg m-2  R/P = ', ratio, '  max |R|/W per step = ', max_rel_step

    if (abs(ratio) > tolerance) then
      write(*,'(A,ES10.2,A,ES10.2)') 'test_water_conservation() --- *** FAILED *** |R/P| = ', &
        abs(ratio), ' exceeds tolerance of ', tolerance
      error stop 'water conservation test failed'
    endif
  end subroutine test_water_conservation


  subroutine test_cloud_fraction_conservation(dt, cldfra, water_tolerance, energy_tolerance, known_issue)
    !! column water and energy budgets with the prognostic cloud fraction scheme
    !!
    !! uses the 59-level MPAS column from test_water_conservation (aerosol aware, with
    !! water vapor added in ascent outside tempo_run so cloud keeps forming) and adds the
    !! cloud fraction inputs: liquid/ice cloud fractions, PBL height, land mask, PBL
    !! subgrid cloud, and a longwave cooling rate. With cldfra = .false. the same column
    !! runs without the scheme as a control.
    !!
    !! water: R = (W_after - W_before) + P per step, as in test_water_conservation, with
    !!   the layer mass from the density TEMPO computes on entry to the step;
    !!   fails if |R/P| > water_tolerance
    !! energy: E = sum(m * (cp*T - Lv(T)*(qc+qr) - Ls*(qi+qs+qg))), with Lv(T) as in TEMPO;
    !!   the residual R_E = (E_after - E_before) - Lv(T_1)*P_liquid - Ls*P_frozen is
    !!   compared with the latent heat of the vapor exchanged, L = sum over steps and
    !!   levels of m*Lv(T)*|dqv|; fails if |R_E/L| > energy_tolerance. E is not exactly
    !!   conserved by TEMPO (moist cp, constant Ls, Lv evaluated at the end of the step),
    !!   so the control run sets the floor for this tolerance
    !!
    !! known_issue = .true. reports a tolerance failure as a known issue instead of
    !! stopping, for cases with an open bug that should still be run and reported
    use module_mp_tempo_params, only : wp, sp, dp, roverrv, rdry, min_qv, cp, lsub, lvap0

    real(wp), intent(in) :: dt
    logical, intent(in) :: cldfra
    real(wp), intent(in) :: water_tolerance, energy_tolerance
    logical, intent(in), optional :: known_issue
    logical :: report_only
    integer :: itimestep
    integer, parameter :: nz = 59
    integer, parameter :: ids = 1, ide = 1, ims = 1, ime = 1, its = 1, ite = 1
    integer, parameter :: jds = 1, jde = 1, jms = 1, jme = 1, jts = 1, jte = 1
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(its:ite, kts:kte, jts:jte) :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr, ng, qb, nc, nwfa, nifa, &
      qcfrac, qifrac, qc_bl, qcfrac_bl, thten_bl, qvten_bl, qcten_bl, qiten_bl, &
      thten_lwrad, thten_swrad
    real(wp), dimension(its:ite, jts:jte) :: hpbl
    integer, dimension(its:ite, jts:jte) :: land
    real(wp), dimension(nz) :: klevs_in, qv_in, qc_in, qr_in, qi_in, qs_in, &
      qg_in, ni_in, nr_in, nc_in, nwfa_in, nifa_in, theta_in, ng_in, volg_in, &
      pressure_in, w_in, dz_in
    real(wp), dimension(nz) :: layer_mass, qv_before
    real(dp) :: water_before, water_after, energy_before, energy_after, &
      precip_liq, precip_frz, precip_sum, water_residual_sum, energy_residual_sum, &
      latent_sum, water_ratio, energy_ratio, cloud_frac_max
    real(wp), parameter :: forcing_gamma = 1.e-6_wp !! vapor added per meter of ascent \([kg\; kg^{-1}\; m^{-1}]\)
    real(wp), parameter :: forcing_w_max = 4._wp !! vertical velocity cap for the forcing \([m\; s^{-1}]\)
    character(len=20) :: dt_string
    character(len=20) :: case_string
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    integer :: io1, io2, k, total_timesteps
    integer, parameter :: integration_time = 1200

    report_only = .false.
    if (present(known_issue)) report_only = known_issue

    write(dt_string, '(I7)') int(dt)
    if (cldfra) then
      case_string = '_cldfra_on'
    else
      case_string = '_cldfra_off'
    endif

    ! read input file
    open(newunit=io1, file='../test/data/mpas_59lev_test.txt', status='old')
    read(io1,*) ! header
    do k = 1, nz
      read(io1,*) klevs_in(k), qv_in(k), qc_in(k), qr_in(k), qi_in(k), qs_in(k), &
          qg_in(k), ni_in(k), nr_in(k), nc_in(k), nwfa_in(k), nifa_in(k), &
          theta_in(k), ng_in(k), volg_in(k), pressure_in(k), w_in(k), dz_in(k)
    end do
    close(io1)

    qv(1,:,1) = qv_in
    t(1,:,1) = theta_in * (pressure_in/100000.)**0.286
    p(1,:,1) = pressure_in
    dz(1,:,1) = dz_in
    w(1,:,1) = w_in
    qc(1,:,1) = qc_in
    nc(1,:,1) = nc_in
    qr(1,:,1) = qr_in
    nr(1,:,1) = nr_in
    qi(1,:,1) = qi_in
    ni(1,:,1) = ni_in
    qs(1,:,1) = qs_in
    qg(1,:,1) = qg_in
    ng(1,:,1) = ng_in
    qb(1,:,1) = volg_in * 1000. ! convert meters^3 -> liters
    nwfa(1,:,1) = nwfa_in
    nifa(1,:,1) = nifa_in

    ! cloud fraction inputs: fractions start at 0 (TEMPO sets them to 1 on the first call),
    ! subgrid PBL cloud in the lowest levels, longwave cooling, and small PBL tendencies
    qcfrac = 0._wp
    qifrac = 0._wp
    qc_bl = 0._wp
    qcfrac_bl = 0._wp
    do k = 2, 6
      qc_bl(1,k,1) = 5.e-5_wp
      qcfrac_bl(1,k,1) = 0.3_wp
    enddo
    thten_bl = 0._wp
    qvten_bl = 1.e-7_wp
    qcten_bl = 1.e-8_wp
    qiten_bl = 0._wp
    thten_lwrad = -2.e-5_wp
    thten_swrad = 0._wp
    hpbl = 1000._wp
    land = 1

    ! set configs
    tempo_cfgs%turn_off_micro_flag = .false.
    tempo_cfgs%semi_sedi_flag = .false.
    tempo_cfgs%cldfra_flag = cldfra

    total_timesteps = int(integration_time/dt)
    precip_sum = 0._dp
    water_residual_sum = 0._dp
    energy_residual_sum = 0._dp
    latent_sum = 0._dp
    cloud_frac_max = 0._dp
    open(newunit=io2, file="cloud_fraction_budget_dt_"//trim(adjustl(dt_string))//trim(case_string)//".txt", &
      status="new", action="write")
    write(io2, '(A)') 'time precip_sum water_residual_sum energy_residual_sum latent_sum'

    do itimestep = 1, total_timesteps
      ! forcing: add water vapor in ascent, outside tempo_run
      do k = 1, nz
        qv(1,k,1) = qv(1,k,1) + forcing_gamma*min(max(w(1,k,1), 0._wp), forcing_w_max)*dt
      enddo

      ! layer mass from the density TEMPO computes on entry to the step
      do k = 1, nz
        layer_mass(k) = roverrv*p(1,k,1) / (rdry*t(1,k,1)*(max(min_qv, qv(1,k,1))+roverrv)) * dz(1,k,1)
      enddo
      qv_before = qv(1,:,1)
      water_before = column_water(layer_mass, qv(1,:,1), qc(1,:,1), qr(1,:,1), &
        qi(1,:,1), qs(1,:,1), qg(1,:,1))
      energy_before = column_energy(layer_mass, t(1,:,1), qc(1,:,1)+qr(1,:,1), &
        qi(1,:,1)+qs(1,:,1)+qg(1,:,1))

      call  tempo_run(tempo_cfgs=tempo_cfgs, itimestep=itimestep, dt=dt, &
                      ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                      jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                      kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                      qc=qc, nc=nc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, ng=ng, qb=qb, &
                      nwfa=nwfa, nifa=nifa, &
                      qc_bl=qc_bl, qcfrac_bl=qcfrac_bl, qcfrac=qcfrac, qifrac=qifrac, hpbl=hpbl, &
                      thten_bl=thten_bl, qvten_bl=qvten_bl, qcten_bl=qcten_bl, qiten_bl=qiten_bl, &
                      thten_lwrad=thten_lwrad, thten_swrad=thten_swrad, land_input=land, &
                      tempo_diags=tempo_driver_diags)

      water_after = column_water(layer_mass, qv(1,:,1), qc(1,:,1), qr(1,:,1), &
        qi(1,:,1), qs(1,:,1), qg(1,:,1))
      energy_after = column_energy(layer_mass, t(1,:,1), qc(1,:,1)+qr(1,:,1), &
        qi(1,:,1)+qs(1,:,1)+qg(1,:,1))

      ! graupel_liquid_equiv_precip includes hail
      precip_liq = tempo_driver_diags%rain_precip(1,1) + tempo_driver_diags%cloud_precip(1,1)
      precip_frz = tempo_driver_diags%ice_liquid_equiv_precip(1,1) + &
        tempo_driver_diags%snow_liquid_equiv_precip(1,1) + tempo_driver_diags%graupel_liquid_equiv_precip(1,1)

      precip_sum = precip_sum + precip_liq + precip_frz
      water_residual_sum = water_residual_sum + (water_after - water_before) + precip_liq + precip_frz
      energy_residual_sum = energy_residual_sum + (energy_after - energy_before) - &
        latent_vap(t(1,1,1))*precip_liq - lsub*precip_frz
      do k = 1, nz
        latent_sum = latent_sum + real(layer_mass(k)*latent_vap(t(1,k,1))*abs(qv(1,k,1) - qv_before(k)), kind=dp)
      enddo
      cloud_frac_max = max(cloud_frac_max, real(maxval(qcfrac), kind=dp), real(maxval(qifrac), kind=dp))
      write(io2, '(I7, 4E16.8)') int(itimestep*dt), precip_sum, water_residual_sum, energy_residual_sum, latent_sum
    enddo
    close(io2)
    tempo_cfgs%cldfra_flag = .false.

    water_ratio = 0._dp
    if (precip_sum > 0._dp) water_ratio = water_residual_sum / precip_sum
    energy_ratio = 0._dp
    if (latent_sum > 0._dp) energy_ratio = energy_residual_sum / latent_sum

    write(*,'(A,I5,A,L1,A,ES11.4,A,ES11.4,A,ES10.2,A,ES11.4,A,ES10.2,A,F6.3)') &
      'test_cloud_fraction_conservation() --- dt = ', int(dt), ' cldfra = ', cldfra, &
      ' P = ', precip_sum, ' kg m-2  R = ', water_residual_sum, ' kg m-2  R/P = ', water_ratio, &
      '  latent = ', latent_sum, ' J m-2  R_E/latent = ', energy_ratio, '  max fraction = ', cloud_frac_max

    if (abs(water_ratio) > water_tolerance) then
      if (report_only) then
        write(*,'(A,ES10.2,A,ES10.2)') 'test_cloud_fraction_conservation() --- *** KNOWN ISSUE *** |R/P| = ', &
          abs(water_ratio), ' exceeds tolerance of ', water_tolerance
      else
        write(*,'(A,ES10.2,A,ES10.2)') 'test_cloud_fraction_conservation() --- *** FAILED *** |R/P| = ', &
          abs(water_ratio), ' exceeds tolerance of ', water_tolerance
        error stop 'cloud fraction water conservation test failed'
      endif
    endif
    if (abs(energy_ratio) > energy_tolerance) then
      if (report_only) then
        write(*,'(A,ES10.2,A,ES10.2)') 'test_cloud_fraction_conservation() --- *** KNOWN ISSUE *** |R_E/latent| = ', &
          abs(energy_ratio), ' exceeds tolerance of ', energy_tolerance
      else
        write(*,'(A,ES10.2,A,ES10.2)') 'test_cloud_fraction_conservation() --- *** FAILED *** |R_E/latent| = ', &
          abs(energy_ratio), ' exceeds tolerance of ', energy_tolerance
        error stop 'cloud fraction energy conservation test failed'
      endif
    endif

  contains

    real(wp) function latent_vap(tk)
      !! latent heat of vaporization as computed in TEMPO
      real(wp), intent(in) :: tk
      latent_vap = lvap0 + (2106._wp - 4218._wp)*(tk - 273.15_wp)
    end function latent_vap

    function column_energy(mass, tk, ql, qf) result(energy)
      !! column-integrated liquid/ice water static energy (without the geopotential term) \([J\; m^{-2}]\)
      real(wp), dimension(:), intent(in) :: mass, tk, ql, qf
      real(dp) :: energy
      integer :: kk
      energy = 0._dp
      do kk = 1, size(mass)
        energy = energy + real(mass(kk)*(cp*tk(kk) - latent_vap(tk(kk))*ql(kk) - lsub*qf(kk)), kind=dp)
      enddo
    end function column_energy
  end subroutine test_cloud_fraction_conservation


  function column_water(layer_mass, qv, qc, qr, qi, qs, qg) result(water)
    !! column-integrated total water \([kg\; m^{-2}]\)
    use module_mp_tempo_params, only : wp, dp

    real(wp), dimension(:), intent(in) :: layer_mass, qv, qc, qr, qi, qs, qg
    real(dp) :: water
    integer :: k

    water = 0._dp
    do k = 1, size(layer_mass)
      water = water + real(layer_mass(k), kind=dp) * &
        real(qv(k) + qc(k) + qr(k) + qi(k) + qs(k) + qg(k), kind=dp)
    enddo
  end function column_water

end module tests