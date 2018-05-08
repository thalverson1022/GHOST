CC=mpif90
FILENAME=ND_HOBasis_v2-1
MODS1=ND_HOBasis_MODS_v2-1
OUTPUT=ghost
TAG=f90
LIBS= -Wl,-rpath,$TACC_MKL_LIB \
    -L$(TACC_SCALAPACK_LIB) -lscalapack \
    $(TACC_SCALAPACK_LIB)/blacs_MPI-LINUX-0.a \
    $(TACC_SCALAPACK_LIB)/blacsF77init_MPI-LINUX-0.a \
    $(TACC_SCALAPACK_LIB)/blacs_MPI-LINUX-0.a \
    -L$(TACC_MKL_LIB) -Wl,--start-group -lmkl_intel_lp64 \
    -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread 
 

all:execute

execute:$(FILENAME).$(TAG) Makefile
	$(CC) $(MODS1).$(TAG) $(FILENAME).$(TAG) -o $(OUTPUT).exe $(LIBS)
clean:
	rm -f $(OUTPUT).exe; 
	rm -f *.mod;

