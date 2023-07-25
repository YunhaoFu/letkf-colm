FC=mpifort
FCFLAGS=-I$(NETCDF)/include
LDFLAGS=-L$(NETCDF)/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -L/opt/intel/oneapi/mkl/2023.1.0/lib/intel64 -lblas -lmkl_rt #-lpthread -ldl # module load mkl/v4.5

objs = precision.o kdtree.o letkf_core.o obs_state.o model_state.o localization.o letkf_main.o letkf.o

letkf.exe : $(objs)
	$(FC) $^ -o $@ $(LDFLAGS)
%.o : %.F90
	$(FC) $(FCFLAGS) -c $^
clean:
	rm -f *.o *.mod letkf.exe
