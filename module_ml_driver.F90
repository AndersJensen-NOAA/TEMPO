module module_ml_driver

  use module_neural_net, only : load_scale_values
  use module_neural_net, only : Dense, init_neural_net, standard_scaler_transform, apply_activation

  implicit none
  integer, parameter :: r8 = selected_real_kind(6)
  integer, parameter :: i8 = selected_int_kind(18)
  
contains

  subroutine ml_driver(input, output)

    real, dimension(:,:), intent(in) :: input
    real, allocatable, dimension(:,:) :: input_transf, nn_input
    integer :: k, n
    integer, parameter :: batch_size = 1
    integer, parameter :: num_outputs = 1
    integer :: iulog
    character(128) :: errstring  ! output status (non-blank for error return)

    type(Dense), allocatable, save :: q_all(:)
    real(r8), dimension(:, :, :), allocatable :: tmp_out, tmp_out_activation
    integer :: i, iout, m, j
    real(r8), dimension(batch_size, num_outputs) :: output_transf
    real(r8), dimension(:), intent(out) :: output
    
    real(r8), allocatable, dimension(:,:), save :: scale_vals
    k = size(input,1) ! model k-dimension (single column)
    n = size(input,2) ! number of ML input variables

    if (.not. allocated(scale_vals)) then
       allocate(scale_vals(n, 2))
    endif
    call load_scale_values("test_scale.txt", n, scale_vals)

    errstring = ''
    iulog = 10
!    write(iulog,*) "Begin loading neural network"
    call init_neural_net("test_nn_real.nc", batch_size, q_all, iulog, errstring)
    if (trim(errstring) /= '') write(iulog,*) 'ERROR in subroutine init_neural_net'

    if (.not. allocated(input_transf)) then
       allocate(input_transf(k, n))
    endif
    call standard_scaler_transform(input, scale_vals, input_transf, errstring)

!    write(*,*) 'transform'
!    write(*,*) input(1,1), input_transf(1,1)
 !   write(*,*) input(1,2), input_transf(1,2)
 !   write(*,*) input(1,3), input_transf(1,3)
    !    write(*,*) input_transf

!    if (.not. allocated(output)) then
!       allocate(output(k))
!    endif
    
    do m = 1, k

       if (.not. allocated(nn_input)) then
          allocate(nn_input(batch_size, n))
       endif
       do j = 1, n
          nn_input(1,j) = input_transf(m,j)
!          write(*,*) m, j, nn_input(1,j)
       enddo
       
       if (.not. allocated(tmp_out)) then
          allocate(tmp_out(size(q_all%output_size), q_all(1)%output_size, num_outputs))
       endif
       if (.not. allocated(tmp_out_activation)) then
          allocate(tmp_out_activation(size(q_all%output_size), q_all(1)%output_size, num_outputs))
       endif
       
       ! Input layer
       tmp_out(1,:,:) = matmul(transpose(q_all(1)%weights), transpose(nn_input)) + &
            reshape(q_all(1)%bias, (/q_all(1)%output_size, num_outputs/))
       call apply_activation(tmp_out(1,:,:), q_all(1)%activation, tmp_out_activation(1,:,:))
       
       ! Hidden layers
       do i = 2, size(q_all%output_size) - 1
          tmp_out(i,:,:) = matmul(transpose(q_all(i)%weights), tmp_out_activation(i-1,:,:)) + &
               reshape(q_all(i)%bias, (/q_all(i)%output_size, num_outputs/))
          call apply_activation(tmp_out(i,:,:), q_all(i)%activation, tmp_out_activation(i,:,:))
       enddo
       
       ! Output layer
       iout = size(q_all%output_size)
       output_transf = matmul(transpose(q_all(iout)%weights), tmp_out_activation(iout-1,:,:)) + &
            reshape(q_all(iout)%bias, (/q_all(iout)%output_size, num_outputs/))
       call apply_activation(output_transf, q_all(iout)%activation, output_transf)

       output(m) = 10**(output_transf(1,1))
       
!        write(*,*) 'OUTPUT', output_transf, 10**(output_transf)/1.e6
    enddo
    
    if (allocated(input_transf)) deallocate(input_transf)
    if (allocated(nn_input)) deallocate(nn_input)
    if (allocated(tmp_out)) deallocate(tmp_out)
    if (allocated(tmp_out_activation)) deallocate(tmp_out_activation)

  end subroutine ml_driver
  
end module module_ml_driver
