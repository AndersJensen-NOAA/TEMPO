FC = gfortran
FCFLAGS = -fdefault-real-8 -fdefault-double-8

SRC=machine.F90 ../module_mp_thompson_params.F90 driver.F90
OBJ=machine.o module_mp_thompson_params.o driver.o

%.o: %.F90
	$(FC) $(FFLAGS) -c $(SRC)

driver: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) 

clean: 
	@rm -f *.mod *.o driver
