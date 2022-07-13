# make [ [all] | [ [WP] [WPx] ] | [clean] | [clean-deps] | [distclean] | [help] ] [EXE=gnu] [AVX=||avx|avx2|avx512] [OMP=||1|]
#   make all : is the same as "make WP WPx" (by default, we will use intel compilers)
#   make WP : no graphics neaded
#   make WPx : graphics (need plplot and X11, set library locations below)
#   make EXE=gnu : means use g++/gfortran instead of intel icpc/ifort, build both WP and WPx (set commands and options below)
#   make WP[x] EXE=gnu : means use g++/gfortran instead of intel icpc/ifort, builds WP[x] (set commands and options below)
#   make clean : removes build subdirectory and executable(s)
#   make clean-deps : removes dependency subdirectory
#   make distclean : remove build directory and executable(s)
#   make help : show options and default

#let's use bash vs. default sh
SHELL=/bin/bash

#prefix for where fftw and plplot installations
LOCAL_DIR0?=.mrm

#assume we want to use intel compilers
#EXE?=intel
#assume we want to use gnu compilers
EXE?=gnu

#which version of avx?
AVX?=

#should fftw be double precsion
FFTD?=

#root build directory
BLD0?=.build

#root depend directory
DEPDIR:=$(BLD0)/dep

#openMP? assume not
OMP?=

ifeq ($(EXE),intel)  #use intel compilers
#intel compilers
	FC=ifort
	CXX=icpc
	CC=icc

#os dependent
ifeq ($(shell uname),Darwin)
          #mac needs these flags for intel wich is the same as -fast minus -static and -xHost
		FFLAGS?=-ipo -mdynamic-no-pic -O3 -no-prec-div -fp-model fast=2
	  #mac puts X11 libraries here
		X11_LIB=-L/usr/X11/lib -lX11
else
         #centos 7 needs these flags for intel which is the same as -fast minus -static and -xHost
		FFLAGS?=-ipo -O3 -no-prec-div -fp-model fast=2
         #centos 7 puts X11 libraries here
		X11_LIB=-L/usr/lib64 -lX11
endif

#prefix where fftw and, optionally, plplot sit
	LOCAL_DIR=$(LOCAL_DIR0)/$(EXE)
	LOCAL_FFTW_DIR=$(LOCAL_DIR)_$(AVX)

#fortran and c++ flags
#	FFLAGS=-fast -axAVX,CORE-AVX2,CORE-AVX512,MIC-AVX512,COMMON-AVX512 -traceback -heap-arrays 2048   #bad idea
	CXXFLAGS?=-fast
	BLD1:=$(EXE)
ifeq ($(AVX),avx)   #for avx
	BLD1:=$(BLD1)_$(AVX)
	override FFLAGS+=-xAVX
	override CXXFLAGS+=-xAVX
else ifeq ($(AVX),avx2)   #for avx2
	BLD1:=$(BLD1)_$(AVX)
	override FFLAGS+=-xCORE-AVX2
	override CXXFLAGS+=-xCORE-AVX2
else ifeq ($(AVX),avx512)   #for avx512
	BLD1:=$(BLD1)_$(AVX)
	override FFLAGS+=-xCORE-AVX512
	override CXXFLAGS+=-xCORE-AVX512
else   #this will optimize to the machine that compiles, not a great idea
	LOCAL_FFTW_DIR=$(LOCAL_DIR)
	override FFLAGS+=-xHost
	override CXXFLAGS+=-xHost
endif

#assume no openMP
ifneq ($(OMP),)
#want to use openMP, for intel this means we need to statically link in libiomp otherwise require library to be on executing host
		override FFLAGS+=-DOMP -qopenmp
		override FFLAGSM=-qopenmp-link=static
		BLD1:=$(BLD1)_omp
endif

ifneq ($(FFTD),)
	BLD1:=$(BLD1)_d
endif
	BLDDIR:=$(BLD0)/$(BLD1)
#	FFLAGS+=-traceback -heap-arrays 2048 -fpp -I$(BLDDIR) -module $(BLDDIR)
	override FFLAGS+=-fpp -module $(BLDDIR)

#older versions of intel need -lstdc++ while newer require -cxxlib
#	FLIBS=-lstdc++
#	FLIBSX=-lstdc++
#  this is ifort (19.0.3.199) - earlier versions  may not support -cxxlib
	override FFLAGSM+=-cxxlib -static-intel

else   #use gnu compilers
	override EXE:=gnu
#gnu compiler
	FC=gfortran
	FCGE430 := $(shell expr $$($(FC) -dumpversion | awk -F. '{for (i=1;i<=NF;i++) s=100*s+$$i; if(s<100) s*=100;if (s<10000) s*=100;print s}') \>= 40300)
ifeq ($(FCGE430),1)
	DEPDIR:=$(BLD0)/dep_gnu
endif

#os dependent
ifeq ($(shell uname),Darwin)
	  #mac needs this version of gcc, if we use it
		CXX=g++-11
		CC=gcc-11
		X11_LIB=-L/usr/X11/lib -lX11
		FLIBS=-lstdc++
		FLIBSX=-lstdc++
else
         #centos 7 needs version 4.9+ of gcc/gfortran to compile plplot, if we use it
		CXX=g++
		CC=gcc
		X11_LIB=-L/usr/lib64 -lX11
		FLIBS=-lstdc++ -static
		FLIBSX=-lstdc++
endif


#prefix where fftw and, optionally, plplot sit
	LOCAL_DIR=$(LOCAL_DIR0)/$(EXE)
	LOCAL_FFTW_DIR=$(LOCAL_DIR)_$(AVX)

#fortran and g++ flags
	FFLAGS?=-O3
	CXXFLAGS?=-O3
	BLD1:=$(EXE)
ifeq ($(AVX),)  #optimize for compiling machine- not a good idea
	LOCAL_FFTW_DIR=$(LOCAL_DIR)
	FFLAGS+=-march=native
	CXXFLAGS+=-march=native
else ifeq ($(AVX),avx512)  #for avx512
	BLD1:=$(BLD1)_$(AVX)
	override FFLAGS+=-m$(AVX)f
	override CXXFLAGS+=-m$(AVX)f
else  #for avx or avx2
	BLD1:=$(BLD1)_$(AVX)
	override FFLAGS+=-m$(AVX)
	override CXXFLAGS+=-m$(AVX)
endif

#assume no openMP
ifneq ($(OMP),)
	override FFLAGS+=-DOMP -fopenmp
	BLD1:=$(BLD1)_omp
endif

ifneq ($(FFTD),)
	BLD1:=$(BLD1)_d
endif
	BLDDIR:=$(BLD0)/$(BLD1)

	#fortran flags
	override FFLAGS+=-ftree-vectorize -fdollar-ok -ffree-line-length-none -fno-range-check -w -cpp -J$(BLDDIR)
	override FFLAGSM=-fsized-deallocation -static-libgcc -static-libgfortran
	#c++ flags
	override CXXFLAGS+=-fpermissive -w -fsized-deallocation 
endif

FFTW_LIBS=-L$(LOCAL_FFTW_DIR)/lib -lfftw3f
FFTWLIB=$(LOCAL_FFTW_DIR)/lib/libfftw3f.a
override FFLAGS+=-I$(BLDDIR) -I$(LOCAL_FFTW_DIR)/include
ifneq ($(FFTD),)
override FFLAGS+=-DFFTD
FFTW_LIBS=-L$(LOCAL_FFTW_DIR)/lib -lfftw3
FFTWLIB=$(LOCAL_FFTW_DIR)/lib/libfftw3.a
endif
PLPLOT_LIB=-L$(LOCAL_DIR)/lib -lplplotfortran -lplplot -lcsirocsa -lqsastime $(X11_LIB) -lpthread
PLPLOT_MOD=-I$(LOCAL_DIR)/lib/fortran/modules/plplot
PLPLOTLIB=$(LOCAL_DIR)/lib/libplplot.a

#everything
SRC:=.
CXXSRC  := $(shell find $(SRC) -name '*.cpp' | sed "s/^\.\///") #get c++ source files
CXXDEPS := $(patsubst %.cpp,$(DEPDIR)/%.d,$(CXXSRC))  #corresponding dependencies files
FSRC    := $(shell find $(SRC) -name '*.f90' | sed "s/^\.\///") #get fortran files except for a few
FDEPS   := $(patsubst %.f90,$(DEPDIR)/%.d,$(FSRC))  #corresponding dependencies files
SRCDIRS := $(shell find $(SRC) \( -name '*.cpp' -o -name '*.f90' \) -exec dirname {} \; | sort -u)  #all directories where source files live

#object files for making WP/WPx
SRC:=src
OBJS := $(patsubst %.cpp,$(BLDDIR)/%.o,$(shell find $(SRC) -name '*.cpp'))  #corresponding object files
OBJS += $(patsubst %.f90,$(BLDDIR)/%.o,$(shell find $(SRC) -name '*.f90' | grep -v grafix.f90))  #corresponding object files

#general compiling commands
COMPILE.f90 = $(FC) $(FFLAGS) 
COMPILE.cxx = $(CXX) $(CXXFLAGS)

.PHONY: all clean distclean help clean-deps blddirs depdirs

#two executables, one w/ and w/o X11/plplot - this is the default target
all : WPx WP

#main executable w/o X11/plplot
WP : $(FFTWLIB) $(CXXDEPS) $(FDEPS) WP_$(BLD1)
	@rm -f $@ && ln -s WP_$(BLD1) $@

#main executable w/ X11/plplot
WPx : $(FFTWLIB) $(PLPLOTLIB) $(CXXDEPS) $(FDEPS) WP_$(BLD1)_X
	@rm -f $@ && ln -s WP_$(BLD1)_X $@

#main executable w/o X11/plplot
WP_$(BLD1) : $(BLDDIR)/src/grafix_nox.o $(OBJS)
	@echo "Creating WP ($@)..."
	$(COMPILE.f90) $(FFLAGSM) -o $@ $(OBJS) $(BLDDIR)/src/grafix_nox.o $(FFTW_LIBS) $(FLIBS)

#main executable w/ X11/plplot
WP_$(BLD1)_X : $(BLDDIR)/src/grafix.o $(OBJS)
	@echo "Creating WPx ($@)..."
	$(COMPILE.f90) $(FFLAGSM) -o $@ $(OBJS) $(BLDDIR)/src/grafix.o $(FFTW_LIBS) $(PLPLOT_LIB) $(FLIBSX)

#compile graphix w/ X11/plplot
$(BLDDIR)/src/grafix_nox.o : src/grafix.f90
	@echo "Compiling $< without plplot/X11 libraries..."
	@$(COMPILE.f90) -c $< -o $@

$(FFTWLIB) :
	@echo "Compiling fftw libraries..."
	./fftw/to_build fftw-3.3.10.tar.gz $(FC) $(CC) $(LOCAL_FFTW_DIR) $(AVX)

$(PLPLOTLIB) :
	@echo "Compiling plplot library..."
	./plplot/to_build 5.15.0 $(FC) $(CC) $(CXX) $(LOCAL_DIR)
	

#compile graphix w/o X11/plplot
$(BLDDIR)/src/grafix.o : src/grafix.f90
	@echo "Compiling $< with plplot/X11 support..."
	@$(COMPILE.f90) -c $(PLPLOT_MOD) -DXWIN $< -o $@

#implicit rule for compiling c++ files
$(BLDDIR)/%.o : %.cpp
$(BLDDIR)/%.o : %.cpp
	@echo "Compiling $<..."
	@$(COMPILE.cxx) -c -o $@ $<

#implicit rule for compiling fortran files
$(BLDDIR)/%.o : %.f90
	@if [ $(notdir $<) != grafix.f90 ]; then \
		echo "Compiling $<..."; \
		$(COMPILE.f90) -c -o $@ $<; \
	fi

#rule to generate dependency files for fortran
$(DEPDIR)/%.d : %.f90 | depdirs blddirs
	@echo "Constructing dependencies for $<..."
ifeq ($(FCGE430),1)
	@./mkdependf90.sh $< '$$(BLDDIR)/' noobj > $@
else
	@./mkdependf90.sh $< '$$(BLDDIR)/' > $@
endif
	@if [ $(notdir $<) == grafix.f90 ]; then \
		rm $@; touch $@; \
	fi
	
#rule to generate dependency files for c++ (note intel has a problem putting an extra $ in dependency file, so we hack it)
$(DEPDIR)/%.d: %.cpp | depdirs blddirs
	@echo "Constructing dependencies for $<..."
	@set -e; rm -f $@; \
	 $(CXX) -MT '$$(BLDDIR)/$(subst .cpp,.o,$<)' -MM $< > $@.T; \
	 sed -e 's,\$$\$$,\$$,' -e 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.T > $@; \
	 rm -f $@.T

#make the build directory- depends on compiler, AVX and OMP
blddirs:
	@[ ! -d "$(BLDDIR)" ] && \
		(echo "Creating build directory $(BLDDIR)..." && \
		(for dir in $(SRCDIRS); do mkdir -p $(BLDDIR)/$$dir > /dev/null; done)) || \
		echo "Reusing build direcotry $(BLDDIR)"

#make the dependency directory- used for all builds
depdirs:
	@[ ! -d "$(DEPDIR)" ] && (echo "Creating dependency directory $(DEPDIR)..." && \
		(for dir in $(SRCDIRS); do mkdir -p $(DEPDIR)/$$dir > /dev/null; done)) || \
		echo "Reusing dependency direcotry $(DEPDIR)"

#only include dependency files if we're actually going to make something (i.e. '', all, WP or WPx)
ifeq ($(MAKECMDGOALS),$(filter $(MAKECMDGOALS),'' WP WPx all))
	@echo "Reading dependencies..."
  -include $(CXXDEPS)
  -include $(FDEPS)
endif

#clean out a particular executable set and build sub-directory
clean :
	@echo "Removing WP_$(BLD1)_X WP_$(BLD1) $(BLDDIR)..."
	@$(RM) -rf WP_$(BLD1)_X WP_$(BLD1) $(BLDDIR)

#clean out dependency files- will be automatically regenerated when needed
clean-deps :
	@echo "Removing $(DEPDIR)..."
	@$(RM) -rf $(DEPDIR)

#wipe out all builds and executables
distclean :
	@echo "Removing WP* $(BLD0)..."
	@$(RM) -rf WP* $(BLD0)
#some help
help :
	@echo "make [all|WP|WPx|help|clean|clean-deps|distclean] [EXE=intel|gnu] [AVX=|avx|avx2|avx512] [OMP=|1] [FFTD=|1]"
	@echo "  default: make all EXE=$(EXE) AVX=$(AVX) OMP=$(OMP)"

