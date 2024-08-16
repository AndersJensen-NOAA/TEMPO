.SUFFIXES: .F90 .o

all: dummy microphysics

dummy:
	echo "****** compiling Thompson-Eidhammer microphysics ******"

OBJS = \
	module_neural_net.o    \
	module_ml_driver.o    \
	module_mp_thompson_params.o    \
	module_mp_thompson_utils.o     \
	module_mp_thompson_main.o      \
	module_mp_thompson.o

microphysics: $(OBJS)
	ar -ru ./../libphys.a $(OBJS)

# DEPENDENCIES:
module_mp_thompson.o: \
	module_mp_thompson_params.o \
	module_mp_thompson_utils.o \
	module_mp_thompson_main.o

module_ml_driver.o: \
	module_neural_net.o

clean:
	$(RM) *.f90 *.o *.mod
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F90.o:
ifeq "$(GEN_F90)" "true"
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -L/apps/netcdf/4.2.1.1-intel/lib -lnetcdff -I/apps/netcdf/4.2.1.1-intel/include -I.. -I../physics_wrf -I../physics_mmm -I../../../framework -I../../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES) -L/apps/netcdf/4.2.1.1-intel/lib -lnetcdff -I/apps/netcdf/4.2.1.1-intel/include -I.. -I../physics_wrf -I../physics_mmm -I../../../framework -I../../../external/esmf_time_f90
endif
