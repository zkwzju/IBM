#### Makefile for channel code #### 
EXEC = chnl.x
TGZFILE = all.tgz
#
SRC = \
modified_windows_fn.F timers_m.F fft_m.F flags_m.F advance.F io.F pstep.F  timers.F viscxyz.F courant.F  rhs.F  \
main.F post_ellip.F initial.F divg.F partial.F  \
stats.F tt_rhs.F post_proc.F nltrms_up.F gen_helmholz.F     \
enlred.F fft2d_new.F idminmax.F    \
 del_fn.F ibm_ellip.F  mt19937.F\
#
INC =  
#
XTRA = Makefile chnl.anz_ini \
       README README_eqns README_gbal README_time README_old
#
OBJS = $(SRC:.F=.o)

########
SRC90= generate_points_m.f90 ellip_common_m.f90 p_dyn_m.f90 global_m.f90 init_m.f90 common_m.f90 precision_m.f90 parser_m.f90  ellipsoid_m.f90 rng_m.f90 rotation_m.f90
OBJ90=$(SRC90:.f90=.o)

########
#
# Note: Requires "module load intel fftw" prior to "make UFHPC=1 MP=1"
# 
ifdef UFHPC
MKLDIR     = $(HPC_MKL_DIR)
MKLROOT    = $(HPC_MKL_DIR)
MKLLIBDIR  = ${MKLDIR}/lib/intel64
MKLINCDIR  = ${MKLDIR}/include
MKLFFTWINC = ${MKLINCDIR}/fftw
LIBS       = -L${MKLLIBDIR} \
             -Wl,--start-group \
	     -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
             -Wl,--end-group 
FLAGS      = -mcmodel=medium -shared-intel -fpp -DIFC -DFFTW3 -I${MKLFFTWINC}
DEBUGFLAGS = -g -C -traceback
OPTFLAGS   = -O2 -mavx -axsse3,sse4.1,sse4.2
MPFLAGS    = -gopenmp -D OPENMP
FCMP       = ifort
FCSP       = ifort
endif
#

ifdef mac
MKLDIR     = /opt/intel/mkl
MKLLIBDIR  = ${MKLDIR}/lib
FFTWDIR    = /usr/local
FFTWLIBDIR = ${FFTWDIR}/lib
LIBS       = -L${FFTWLIBDIR} -lfftw3 -L${MKLLIBDIR} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
FLAGS      =  -shared-intel -fpp -traceback -DIFC -DFFTW3 -I${FFTWDIR}/include
DEBUGFLAGS = -g -C -traceback
OPTFLAGS   = -O2
MPFLAGS    = -qopenmp -D OPENMP
FCMP       = ifort
FCSP       = ifort
endif

TARFILE    = chnl.tar
#
#ifdef MP
  FC = $(FCMP)
  FLAGS += $(MPFLAGS)
#else
#  FC = $(FCSP)
#endif
ifdef DEBUG
  FLAGS += $(DEBUGFLAGS)
else
  FLAGS += $(OPTFLAGS)
endif

%.o: %.F $(INC)
	$(FC) -c $(FLAGS) $<


%.o: %.f90 $(INC)
	$(FC) -c $(FLAGS) $<


$(EXEC): $(OBJ90) $(OBJS)
	$(FC) -o $@ $(FLAGS) $^ $(LIBS)

dist:
	-rm -f $(TARFILE).gz
	tar -cf $(TARFILE) $(SRC) $(INC) $(XTRA)
	gzip $(TARFILE)

get:
	scp mcantero@cee-zze6.cee.uiuc.edu:~/research/spectral_code/front_version/code/code_mixed9.0_aniso/$(TGZFILE) .
	tar -zxvf $(TGZFILE)
	rm -f $(TGZFILE)

clean:
	-rm -f $(OBJS) $(OBJ90) *.mod


modified_windows_fn.o : common_m.o modified_windows_fn.F
ellip_common_m.o : common_m.o rng_m.o
generate_points_m.o : ellip_common_m.o 
ellipsoid_m.o : generate_points_m.o common_m.o parser_m.o rng_m.o ellipsoid_m.f90
initial.o     : ellipsoid_m.o initial.F
rng_m.o	      : rng_m.f90 mt19937.o 
del_fn.o      : ellipsoid_m.o del_fn.F
parser_m.o    : parser_m.f90 precision_m.o
rho.F	      : ibm_ellip.o rhs.F
ibm_ellip.o   :	rotation_m.o ellipsoid_m.o ibm_ellip.F
post_ellip.o  : ellipsoid_m.o
global_m.o    : parser_m.o
common_m.o    : global_m.o flags_m.o fft_m.o timers_m.o fft_m.o common_m.f90
fft_m.o	      : global_m.o fft_m.F
init_m.o      : common_m.o init_m.f90
p_dyn_m.o     : rotation_m.o common_m.o ellip_common_m.o p_dyn_m.f90
advance.o     : common_m.o ellip_common_m.o p_dyn_m.o
main.o	      : init_m.o ellipsoid_m.o main.F

