.SUFFIXES: .o .F


F90  = /usr/local/mpich2-1.0.7/bin/mpif90
## F90 = ifort

# FITSIO Linking
FITS_LINK = -L/shared/cfitsio/lib -lcfitsio


# MKL Linking Set up as per the users guide (see users guide chapters 5 and 9) N.F. 2/17/10

MKLPATH = /central/intel/mkl_10.1.0.015
MKL_START = -L${MKLPATH}/lib/em64t -I${MKLPATH}/include -Wl,--start-group
MKL_END = -static_mpi -Wl,--end-group -lpthread -lm
MKL_END = -Wl,--end-group -lpthread -lm
CLUSTER =  ${MKLPATH}/lib/em64t/libmkl_scalapack_lp64.a
BLACS = ${MKLPATH}/lib/em64t/libmkl_blacs_intelmpi_lp64.a
CORE_LAPACK = ${MKLPATH}/lib/em64t/libmkl_intel_lp64.a
CLUSTERI =  ${MKLPATH}/lib/em64t/libmkl_scalapack_ilp64.a
BLACSI = ${MKLPATH}/lib/em64t/libmkl_blacs_intelmpi_ilp64.a
CORE_LAPACKI = ${MKLPATH}/lib/em64t/libmkl_intel_ilp64.a
CORE_KERN1 = ${MKLPATH}/lib/em64t/libmkl_sequential.a
CORE_KERN2 = ${MKLPATH}/lib/em64t/libmkl_core.a
MKL_LINK = $(MKL_START) $(CLUSTER) $(BLACS) $(CORE_LAPACK) $(CORE_KERN1) $(CORE_KERN2) $(MKL_END)
MKL_LINKI = $(MKL_START) $(CLUSTERI) $(BLACSI) $(CORE_LAPACKI) $(CORE_KERN1) $(CORE_KERN2) $(MKL_END)

# Add it all together
LIBFLAGS = $(MKL_LINK) $(FITS_LINK)
## LIBFLAGS = $(MKL_LINKI) $(FITS_LINK)
## use the commented flags for debugging
F90FLAGS = -g -traceback -CB -convert big_endian -fpp -FR -m64 -xO -msse3 -O3 -r8 -I../include 
# one of the following flags is absolutely necessary, but I haven't figured out which one
F90FLAGS = -convert big_endian -fpp -FR -m64 -xO -msse3 -O3 -ip -ipo -r8 -I../include
#F90FLAGS = -convert big_endian -fpp -FR -I../include
F90FLAGS_OMP = -openmp -convert big_endian -fpp -FR -m64 -xO -msse3 -O0 -r8 -I../include 





arrdi : Main.o Run_Info.o Matrix_Magic.o Fitslib_NF.o Solution_Grid.o Region_Info.o Data_Management.o Kernel_Management.o IO_Routines.o
	$(F90) $(F90FLAGS) *.o -o $@ $(LIBFLAGS)

.F.o :
	$(F90) $(F90FLAGS) -c $<

clean :
	rm -f *.o *.mod 

Main.o: Fitslib_NF.o Run_Info.o Matrix_Magic.o Solution_Grid.o Data_Management.o Region_Info.o
Run_Info.o:  
Matrix_Magic.o:  Run_Info.o Fitslib_NF.o Data_Management.o Region_Info.o Solution_Grid.o Kernel_Management.o
Solution_Grid.o:  Fitslib_NF.o Kernel_Management.o Run_Info.o
Region_Info.o: Solution_Grid.o
Data_Management.o: Solution_Grid.o Region_Info.o IO_Routines.o Kernel_Management.o Run_Info.o
Kernel_Management.o: IO_Routines.o Fitslib_NF.o
IO_Routines.o:
