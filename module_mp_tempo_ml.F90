module module_mp_tempo_ml

  implicit none
  private
  public :: tempo_create_nn, predict_nc, ml_data
    
  type ml_data
     integer :: input_size
     integer :: output_size
     integer :: node_size
     real, allocatable, dimension(:) :: transformMean
     real, allocatable, dimension(:) :: transformVar
     real, allocatable, dimension(:,:) :: weights00
     real, allocatable, dimension(:) :: bias00
     real, allocatable, dimension(:,:) :: weights01
     real, allocatable, dimension(:) :: bias01
  end type ml_data

contains

  subroutine tempo_create_nn(tempo_ml_data_in, tempo_ml_data_out)

    logical, save :: nnNotInitialized = .true.
    type(ml_data), intent(in), optional :: tempo_ml_data_in
    type(ml_data), intent(out), optional :: tempo_ml_data_out
    type(ml_data), save :: tempo_ml_data_save

    if (nnNotInitialized) then
       if (present(tempo_ml_data_in)) then
          tempo_ml_data_save = tempo_ml_data_in
       endif
       nnNotInitialized = .false.
    endif

    if (present(tempo_ml_data_out)) then
       tempo_ml_data_out = tempo_ml_data_save
    endif
    
  end subroutine tempo_create_nn

!---------------------------------------

  real function predict_nc(qc, qr, qi, qs, pres, temp, w)

    real, intent(in) :: qc, qr, qi, qs, pres, temp, w
    type(ml_data) :: get_ml_data
    integer, parameter :: input_rows = 1
    real, allocatable, dimension(:,:) :: NNinput, NNinput_transformed
    real, allocatable, dimension(:,:) :: output00, output00Activ
    real, allocatable, dimension(:,:) :: output01, output01Activ
    real, parameter :: logMin = 4.
    real, parameter :: logMax = 9.3010299957
    real :: predictExp
    
    ! Get NN data
    call tempo_create_nn(tempo_ml_data_out=get_ml_data)

    if (.not. allocated(NNinput)) then
       allocate(NNinput(input_rows, get_ml_data%input_size))
    endif
    if (.not. allocated(NNinput_transformed)) then
       allocate(NNinput_transformed(input_rows, get_ml_data%input_size))
    endif
    
    ! Collect input data
    NNinput(1,1) = qc
    NNinput(1,2) = qr
    NNinput(1,3) = qi
    NNinput(1,4) = qs
    NNinput(1,5) = pres
    NNinput(1,6) = temp
    NNinput(1,7) = w

    ! Transform input data
    call standard_scaler_transform(mean=get_ml_data%transformMean, var=get_ml_data%transformVar, &
         raw_data=NNinput, transformed_data=NNinput_transformed)
      
    if (.not. allocated(output00)) then
       allocate(output00(get_ml_data%node_size, get_ml_data%output_size))
    endif
    if (.not. allocated(output00Activ)) then
       allocate(output00Activ(get_ml_data%node_size, get_ml_data%output_size))
    endif
       
    output00 = matmul(transpose(get_ml_data%weights00), transpose(NNinput_transformed)) + &
         reshape(get_ml_data%bias00, (/size(get_ml_data%bias00), get_ml_data%output_size/))

    call relu_activation(input=output00, output=output00Activ)
    
    if (.not. allocated(output01)) then
       allocate(output01(get_ml_data%output_size, get_ml_data%output_size))
    endif
    if (.not. allocated(output01Activ)) then
       allocate(output01Activ(get_ml_data%output_size, get_ml_data%output_size))
    endif
    
    output01 = matmul(transpose(get_ml_data%weights01), output00Activ) + &
         reshape(get_ml_data%bias01, (/size(get_ml_data%bias01), get_ml_data%output_size/))
    
    call relu_activation(input=output01, output=output01Activ)
    
    predictExp = min(logMax, max(logMin, output01Activ(1,1)))
    predict_nc = 10.**predictExp
    
  end function predict_nc


  subroutine standard_scaler_transform(mean, var, raw_data, transformed_data)
    real, intent(in) :: raw_data(:,:)
    real, intent(in) :: mean(:), var(:)
    real, intent(out) :: transformed_data(size(raw_data, 1), size(raw_data, 2))
    integer :: i
    
    do i = 1, size(raw_data, 2)
       transformed_data(:,i) = (raw_data(:,i) - mean(i)) / sqrt(var(i))
    end do
  end subroutine standard_scaler_transform

  
  subroutine relu_activation(input, output)

    real, dimension(:, :), intent(in) :: input
    real, dimension(size(input, 1), size(input, 2)), intent(out) :: output
    integer :: i, j
    
    do i = 1, size(input, 1)
       do j = 1, size(input,2)
          output(i, j) = max(input(i,j), 0.)
       end do
    end do

  end subroutine relu_activation
  
end module module_mp_tempo_ml
