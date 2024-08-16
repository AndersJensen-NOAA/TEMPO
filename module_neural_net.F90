module module_neural_net
    use netcdf

    implicit none
    integer, parameter :: r8 = selected_real_kind(6)
    type Dense
        integer :: input_size
        integer :: output_size
        integer :: batch_size
        integer :: activation
        real(kind=r8), allocatable :: weights(:, :)
        real(kind=r8), allocatable :: bias(:)
    end type Dense

    type DenseData
        real(kind=r8), allocatable :: input(:, :)
        real(kind=r8), allocatable :: output(:, :)
    end type DenseData

contains

    subroutine apply_activation(input, activation_type, output)
        ! Description: Apply a nonlinear activation function to a given array of input values.
        !
        ! Inputs:
        ! input: A 2D array
        ! activation_type: string describing which activation is being applied. If the activation
        !       type does not match any of the available options, the linear activation is applied.
        !       Currently supported activations are:
        !           relu
        !           elu
        !           selu
        !           sigmoid
        !           tanh
        !           softmax
        !           linear
        ! Output:
        ! output: Array of the same dimensions as input with the nonlinear activation applied.
        real(kind=r8), dimension(:, :), intent(in) :: input
        integer, intent(in) :: activation_type
        real(kind=r8), dimension(size(input, 1), size(input, 2)), intent(out) :: output

        real(kind=r8), dimension(size(input, 1)) :: softmax_sum
        real(kind=r8), parameter :: selu_alpha = 1.6732_r8
        real(kind=r8), parameter :: selu_lambda = 1.0507_r8
        real(kind=r8), parameter :: zero = 0.0_r8
        integer :: i, j
        select case (activation_type)
            case (0)
                output = input
            case (1)
                do i=1,size(input, 1)
                    do j=1, size(input,2)
                        output(i, j) = dmax1(DBLE(input(i, j)), zero)
                    end do
                end do
            case default
                output = input
        end select
    end subroutine apply_activation

    subroutine init_neural_net(filename, batch_size, neural_net_model, iulog, errstring)
        ! init_neuralnet
        ! Description: Loads dense neural network weights from a netCDF file and builds an array of
        ! Dense types from the weights and activations.
        !
        ! Input:
        ! filename: Full path to the netCDF file
        ! batch_size: number of items in single batch. Used to set intermediate array sizes.
        !
        ! Output:
        ! neural_net_model (output): array of Dense layers composing a densely connected neural network
        !
        character(len=*), intent(in) :: filename
        integer, intent(in) :: batch_size
        type(Dense), allocatable, intent(out) :: neural_net_model(:)
        integer,          intent(in)  :: iulog
        character(128),   intent(out) :: errstring  ! output status (non-blank for error return)

        integer :: ncid, num_layers_id, num_layers
        integer :: layer_names_var_id, i, layer_in_dimid, layer_out_dimid
        integer :: layer_in_dim, layer_out_dim
        integer :: layer_weight_var_id
        integer :: layer_bias_var_id

        character (len=8), allocatable :: layer_names(:)
        character (len=10) :: num_layers_dim_name = "num_layers"
        character (len=11) :: layer_name_var = "layer_names"
        character (len=11) :: layer_in_dim_name
        character (len=12) :: layer_out_dim_name
        character (len=10) :: activation_name
        real (kind=r8), allocatable :: temp_weights(:, :)

        errstring = ''
        ! Open netCDF file
        call check(nf90_open(filename, nf90_nowrite, ncid),errstring)
        if (trim(errstring) /= '') return
        ! Get the number of layers in the neural network
        call check(nf90_inq_dimid(ncid, num_layers_dim_name, num_layers_id),errstring)
        if (trim(errstring) /= '') return
        call check(nf90_inquire_dimension(ncid, num_layers_id, &
                                          num_layers_dim_name, num_layers),errstring)
        if (trim(errstring) /= '') return
        call check(nf90_inq_varid(ncid, layer_name_var, layer_names_var_id),errstring)
        if (trim(errstring) /= '') return
        allocate(layer_names(num_layers))
        call check(nf90_get_var(ncid, layer_names_var_id, layer_names),errstring)
        if (trim(errstring) /= '') return
        write(iulog,*) "load neural network " // filename
        allocate(neural_net_model(1:num_layers))
        ! Loop through each layer and load the weights, bias term, and activation function
        do i=1, num_layers
            layer_in_dim_name = trim(layer_names(i)) // "_in"
            layer_out_dim_name = trim(layer_names(i)) // "_out"
            layer_in_dimid = -1
            ! Get layer input and output dimensions
            call check(nf90_inq_dimid(ncid, trim(layer_in_dim_name), layer_in_dimid),errstring)
            if (trim(errstring) /= '') return
            call check(nf90_inquire_dimension(ncid, layer_in_dimid, layer_in_dim_name, layer_in_dim),errstring)
            if (trim(errstring) /= '') return
            call check(nf90_inq_dimid(ncid, trim(layer_out_dim_name), layer_out_dimid),errstring)
            if (trim(errstring) /= '') return
            call check(nf90_inquire_dimension(ncid, layer_out_dimid, layer_out_dim_name, layer_out_dim),errstring)
            if (trim(errstring) /= '') return
            call check(nf90_inq_varid(ncid, trim(layer_names(i)) // "_weights", &
                                      layer_weight_var_id),errstring)
            if (trim(errstring) /= '') return
            call check(nf90_inq_varid(ncid, trim(layer_names(i)) // "_bias", &
                                      layer_bias_var_id),errstring)
            if (trim(errstring) /= '') return
            neural_net_model(i)%input_size = layer_in_dim
            neural_net_model(i)%output_size = layer_out_dim
            neural_net_model(i)%batch_size = batch_size
            ! Fortran loads 2D arrays in the opposite order from Python/C, so I
            ! first load the data into a temporary array and then apply the
            ! transpose operation to copy the weights into the Dense layer
            allocate(neural_net_model(i)%weights(layer_in_dim, layer_out_dim))
            allocate(temp_weights(layer_out_dim, layer_in_dim))

            call check(nf90_get_var(ncid, layer_weight_var_id, &
                                    temp_weights),errstring)
            if (trim(errstring) /= '') return
            neural_net_model(i)%weights = transpose(temp_weights)
            deallocate(temp_weights)
            ! Load the bias weights
            allocate(neural_net_model(i)%bias(layer_out_dim))
            call check(nf90_get_var(ncid, layer_bias_var_id, &
                                    neural_net_model(i)%bias),errstring)
            if (trim(errstring) /= '') return
            ! Get the name of the activation function, which is stored as an attribute of the weights variable
            call check(nf90_get_att(ncid, layer_weight_var_id, "activation", &
                                    activation_name),errstring)
            if (trim(errstring) /= '') return
            select case (trim(activation_name))
                case ("linear")
                    neural_net_model(i)%activation = 0
                case ("relu")
                    neural_net_model(i)%activation = 1
                case ("sigmoid")
                    neural_net_model(i)%activation = 2
                case ("elu")
                    neural_net_model(i)%activation = 3
                case ("selu")
                    neural_net_model(i)%activation = 4
                case ("tanh")
                    neural_net_model(i)%activation = 5
                case ("softmax")
                    neural_net_model(i)%activation = 6
                case default
                    neural_net_model(i)%activation = 7
            end select
        end do
        call check(nf90_close(ncid),errstring)
        if (trim(errstring) /= '') return

    end subroutine init_neural_net

    subroutine load_quantile_scale_values(filename, scale_values, iulog, errstring)
        character(len = *), intent(in) :: filename
        real(kind = r8), allocatable, intent(out) :: scale_values(:, :)
        integer,          intent(in)  :: iulog
        character(128),   intent(out) :: errstring  ! output status (non-blank for error return)

        real(kind = r8), allocatable :: temp_scale_values(:, :)
        character(len=8) :: quantile_dim_name = "quantile"
        character(len=7) :: column_dim_name = "column"
        character(len=9) :: ref_var_name = "reference"
        character(len=9) :: quant_var_name = "quantiles"
        integer :: ncid, quantile_id, column_id, quantile_dim, column_dim, ref_var_id, quant_var_id
        
        errstring = ''

        call check(nf90_open(filename, nf90_nowrite, ncid),errstring)
        if (trim(errstring) /= '') return
        call check(nf90_inq_dimid(ncid, quantile_dim_name, quantile_id),errstring)
        if (trim(errstring) /= '') return
        call check(nf90_inq_dimid(ncid, column_dim_name, column_id),errstring)
        if (trim(errstring) /= '') return
        call check(nf90_inquire_dimension(ncid, quantile_id, &
                quantile_dim_name, quantile_dim),errstring)
        if (trim(errstring) /= '') return
        call check(nf90_inquire_dimension(ncid, column_id, &
                column_dim_name, column_dim),errstring)
        if (trim(errstring) /= '') return
        allocate(scale_values(quantile_dim, column_dim + 1))
        allocate(temp_scale_values(column_dim + 1, quantile_dim))
        call check(nf90_inq_varid(ncid, ref_var_name, ref_var_id),errstring)
        if (trim(errstring) /= '') return
        write(iulog,*) "load ref var"
        call check(nf90_get_var(ncid, ref_var_id, temp_scale_values(1, :)),errstring)
        if (trim(errstring) /= '') return
        call check(nf90_inq_varid(ncid, quant_var_name, quant_var_id),errstring)
        if (trim(errstring) /= '') return
        write(iulog,*) "load quant var"
        call check(nf90_get_var(ncid, quant_var_id, temp_scale_values(2:column_dim + 1, :)),errstring)
        if (trim(errstring) /= '') return
        scale_values = transpose(temp_scale_values)
        call check(nf90_close(ncid),errstring)
        if (trim(errstring) /= '') return
    end subroutine load_quantile_scale_values

    subroutine linear_interp(x_in, xs, ys, y_in)
        real(kind = r8), dimension(:), intent(in) :: x_in
        real(kind = r8), dimension(:), intent(in) :: xs
        real(kind = r8), dimension(:), intent(in) :: ys
        real(kind = r8), dimension(size(x_in, 1)), intent(out) :: y_in
        integer :: i, j, x_in_size, xs_size, x_pos
        x_in_size = size(x_in, 1)
        xs_size = size(xs, 1)
        do i = 1, x_in_size
            if (x_in(i) <= xs(1)) then
                y_in(i) = ys(1)
            else if (x_in(i) >= xs(xs_size)) then
                y_in(i) = ys(xs_size)
            else
                j = 1
                do while (xs(j) < x_in(i))
                    j = j + 1
                end do
                y_in(i) = (ys(j - 1) * (xs(j) - x_in(i)) + ys(j) * (x_in(i) - xs(j - 1))) / (xs(j) - xs(j - 1))
            end if
        end do
    end subroutine linear_interp

    subroutine quantile_transform(x_inputs, scale_values, x_transformed)
        real(kind = r8), dimension(:, :), intent(in) :: x_inputs
        real(kind = r8), dimension(:, :), intent(in) :: scale_values
        real(kind = r8), dimension(size(x_inputs, 1), size(x_inputs, 2)), intent(out) :: x_transformed
        integer :: j, x_size, scale_size
        x_size = size(x_inputs, 1)
        scale_size = size(scale_values, 1)
        do j = 1, size(x_inputs, 2)
            call linear_interp(x_inputs(:, j), scale_values(:, j + 1), &
                    scale_values(:, 1), x_transformed(:, j))
        end do
    end subroutine quantile_transform

    subroutine quantile_inv_transform(x_inputs, scale_values, x_transformed)
        real(kind = r8), dimension(:, :), intent(in) :: x_inputs
        real(kind = r8), dimension(:, :), intent(in) :: scale_values
        real(kind = r8), dimension(size(x_inputs, 1), size(x_inputs, 2)), intent(out) :: x_transformed
        integer :: j, x_size, scale_size
        x_size = size(x_inputs, 1)
        scale_size = size(scale_values, 1)
        do j = 1, size(x_inputs, 2)
            call linear_interp(x_inputs(:, j), scale_values(:, 1), scale_values(:, j + 1), x_transformed(:, j))
        end do
    end subroutine quantile_inv_transform

    subroutine standard_scaler_transform(input_data, scale_values, transformed_data, errstring)
        ! Perform z-score normalization of input_data table. Equivalent to scikit-learn StandardScaler.
        !
        ! Inputs:
        !   input_data: 2D array where rows are examples and columns are variables
        !   scale_values: 2D array where rows are the input variables and columns are mean and standard deviation
        ! Output:
        !   transformed_data: 2D array with the same shape as input_data containing the transformed values.
        real(r8), intent(in) :: input_data(:, :)
        real(r8), intent(in) :: scale_values(:, :)
        real(r8), intent(out) :: transformed_data(size(input_data, 1), size(input_data, 2))
        character(128),   intent(out) :: errstring  ! output status (non-blank for error return)
        integer :: i
 
        errstring = ''
        if (size(input_data, 2) /= size(scale_values, 1)) then
            write(errstring,*) "Size mismatch between input data and scale values", size(input_data, 2), size(scale_values, 1)
            return
        end if
        do i=1, size(input_data, 2)
            transformed_data(:, i) = (input_data(:, i) - scale_values(i, 1)) / sqrt(scale_values(i, 2))
        end do
    end subroutine standard_scaler_transform

    subroutine load_scale_values(filename, num_inputs, scale_values)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: num_inputs
        real(r8), intent(out) :: scale_values(num_inputs, 2)
        character(len=40) :: row_name
        integer :: isu, i
        isu = 2
        open(isu, file=filename, access="sequential", form="formatted")
        read(isu, "(A)")
        do i=1, num_inputs
            read(isu, *) scale_values(i, 1), scale_values(i, 2)
        end do
        close(isu)
    end subroutine load_scale_values


    subroutine standard_scaler_inverse_transform(input_data, scale_values, transformed_data, errstring)
        ! Perform inverse z-score normalization of input_data table. Equivalent to scikit-learn StandardScaler.
        !
        ! Inputs:
        !   input_data: 2D array where rows are examples and columns are variables
        !   scale_values: 2D array where rows are the input variables and columns are mean and standard deviation
        ! Output:
        !   transformed_data: 2D array with the same shape as input_data containing the transformed values.
        real(r8), intent(in) :: input_data(:, :)
        real(r8), intent(in) :: scale_values(:, :)
        real(r8), intent(out) :: transformed_data(size(input_data, 1), size(input_data, 2))
        character(128),   intent(out) :: errstring  ! output status (non-blank for error return)
        integer :: i
        if (size(input_data, 2) /= size(scale_values, 1)) then
            write(errstring,*) "Size mismatch between input data and scale values", size(input_data, 2), size(scale_values, 1)
            return
        end if
        do i=1, size(input_data, 2)
            transformed_data(:, i) = input_data(:, i) * scale_values(i, 2) + scale_values(i, 1)
        end do
    end subroutine standard_scaler_inverse_transform

    subroutine minmax_scaler_transform(input_data, scale_values, transformed_data, errstring)
        ! Perform min-max scaling of input_data table. Equivalent to scikit-learn MinMaxScaler.
        !
        ! Inputs:
        !   input_data: 2D array where rows are examples and columns are variables
        !   scale_values: 2D array where rows are the input variables and columns are min and max.
        ! Output:
        !   transformed_data: 2D array with the same shape as input_data containing the transformed values.
        real(r8), intent(in) :: input_data(:, :)
        real(r8), intent(in) :: scale_values(:, :)
        real(r8), intent(out) :: transformed_data(size(input_data, 1), size(input_data, 2))
        character(128),   intent(out) :: errstring  ! output status (non-blank for error return)

        integer :: i
        if (size(input_data, 2) /= size(scale_values, 1)) then
            write(errstring,*) "Size mismatch between input data and scale values", size(input_data, 2), size(scale_values, 1)
            return
        end if
        do i=1, size(input_data, 2)
            transformed_data(:, i) = (input_data(:, i) - scale_values(i, 1)) / (scale_values(i, 2) - scale_values(i ,1))
        end do
    end subroutine minmax_scaler_transform

    subroutine minmax_scaler_inverse_transform(input_data, scale_values, transformed_data, errstring)
        ! Perform inverse min-max scaling of input_data table. Equivalent to scikit-learn MinMaxScaler.
        !
        ! Inputs:
        !   input_data: 2D array where rows are examples and columns are variables
        !   scale_values: 2D array where rows are the input variables and columns are min and max.
        ! Output:
        !   transformed_data: 2D array with the same shape as input_data containing the transformed values.
        real(r8), intent(in) :: input_data(:, :)
        real(r8), intent(in) :: scale_values(:, :)
        real(r8), intent(out) :: transformed_data(size(input_data, 1), size(input_data, 2))
        character(128),   intent(out) :: errstring  ! output status (non-blank for error return)

        integer :: i
        if (size(input_data, 2) /= size(scale_values, 1)) then
            write(errstring,*) "Size mismatch between input data and scale values", size(input_data, 2), size(scale_values, 1)
            return
        end if
        do i=1, size(input_data, 2)
            transformed_data(:, i) = input_data(:, i) * (scale_values(i, 2) - scale_values(i ,1)) + scale_values(i, 1)
        end do
    end subroutine minmax_scaler_inverse_transform

    subroutine check(status, errstring)
        ! Check for netCDF errors
        integer, intent ( in) :: status
        character(128),   intent(out) :: errstring  ! output status (non-blank for error return)

        errstring = ''
        if(status /= nf90_noerr) then
          errstring = trim(nf90_strerror(status))
        end if
    end subroutine check

end module module_neural_net
