ifeq ($(FCOMP),gfortran)
FCFLAGS = -O2 -cpp -freal-4-real-8 -march=native -flto -mavx
FCDEBUG = -g -fbacktrace -fcheck=all -fbounds-check -ffpe-trap=invalid,overflow,underflow,denormal
FCBUILD = -Wall -Wextra -pedantic -std=f2008
FCOPENMP = -fopenmp
else ifeq ($(FCOMP),mpifort)
FCFLAGS = -O2 -cpp -freal-4-real-8 -march=native -flto -mavx -DMPI
FCDEBUG = -g -fbacktrace -fcheck=all -fbounds-check -ffpe-trap=invalid,overflow,underflow,denormal
FCBUILD = -Wall -Wextra -pedantic -std=f2008
else ifeq ($(FCOMP),mpiifort)
FCFLAGS = -O2 -cpp -Dintel -r8 -xHOST -Dintel -ipo -assume byterecl
FCDEBUG = -g -debug -traceback -check all -CB -CU -CA -fpe0 -debug
FCBUILD = -warn all -stand f08
endif
PROGRAM =  mcgrid

SRCS =      random_mod.f90 \
		    utils.f90 \
            constants.f90 \
            stackMod.f90 \
            vector_class.f90 \
            sdfsMod.f90 \
            ch_opt.f90 \
            surfaces.f90 \
            grid.f90 \
            iarray.f90 \
            subs.f90 \
		    photon.f90 \
            inttau2.f90 \
            stokes.f90 \
            writer.f90 \
            main.f90 
            
OBJECTS = $(SRCS:.f90=.o)

all:	$(PROGRAM)
mp:		FCFLAGS += $(FCOPENMP)
mp:		$(PROGRAM)
# mpdebg:	FCFLAGS += $(FCDEBUG)
# mpdebg:	FCFLAGS += $(FCOPENMP)
# mpdebg:	$(PROGRAM)
debug:	FCFLAGS += $(FCDEBUG)
debug:	$(PROGRAM)
build:	FCFLAGS += $(FCBUILD)
build:	$(PROGRAM)

$(PROGRAM):	$(OBJECTS)
	$(FCOMP) $(FCFLAGS)	-o $@ $^ 
	
%.o:  %.f90
	$(FCOMP)  $(FCFLAGS) -c $<

.PHONY:	clean

clean:
	rm -f *.o *.mod *.MOD mcgrid
