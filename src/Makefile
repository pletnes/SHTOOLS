######################################################################
#
#	make all
#		Compile program in the current directory. Optionally, one
#		can specify the parameters F95="my compiler" and 
#		F95FLAGS="my compiler flags". The default is to use "f95".
#
#	make all2
#		Compile program in the current directory. Optionally, one
#		can specify the parameters F95="my compiler" and 
#		F95FLAGS="my compiler flags". LAPACK subroutine names have
#		an underscore appended to them in order to use FFTW and LAPACK
#		libraries with conflicting underscore convention. The default 
#		is to use "f95".
#
#	make clean
#		Remove the lib and module files in the LOCAL directory.
#
#
#	Written by Mark Wieczorek (July 2006).
#
#######################################################################

#############################################
# This block generates the python interface.
# Added by Paul Anton Letnes, May 2011
# This is a work in progress!
#############################################

python: pyshtools.so

# Location of fftw3
LIB_LOCATION = /Users/paulanto/usr/local/lib
LIBS = -llapack -lfftw3
GFORTRAN_FLAGS = -fPIC -fbounds-check -ftrapv -fno-automatic -Wall -Wuninitialized -O0 -pedantic -fimplicit-none -fno-second-underscore
pyshtools.so: pyshtools.f95 $(PROG)
	gfortran -c SHTOOLS.f95
	gfortran -c $(GFORTRAN_FLAGS) pyshtools.f95
	f2py -c -m pyshtools pyshtools.f95 shtools.o $(PROG) -L$(LIB_LOCATION) $(LIBS) --fcompiler=gnu95 -DF2PY_REPORT_ON_ARRAY_COPY=1 --f90flags='$(GFORTRAN_FLAGS)'

shtools.so: shtools.pyf *.f95
	f2py -c shtools.pyf PlmBar.f95 SHTOOLS.f95

#############################################

SHELL=/bin/tcsh

#F95 = f95
F95 = gfortran

ifeq ($(F95),f95)
	# Default Absoft Pro Fortran flags
	#F95FLAGS ?= -O3 -YEXT_NAMES=LCS
	#F95FLAGS ?= -O3 -YEXT_NAMES=LCS -YEXT_SFX=_
	F95FLAGS ?= -O3 -YEXT_NAMES=LCS -YEXT_SFX=_ -m64
endif

ifeq ($(F95),g95)
	# Default g95 flags
	#F95FLAGS ?= -O3 -fno-underscoring
	F95FLAGS ?= -O3 -fno-second-underscore
	#F95FLAGS ?= -O3 -fno-second-underscore -m64
endif


LIBTOOL = libtool
LIBTOOLFLAGS = -static
AR = ar
ARFLAGS = -r
RLIB = ranlib
RLIBFLAGS =
LDFLAGS = -s
MAKE = make


PROG =	../lib/libSHTOOLS2.5.a

SRCS0 =	CilmPlus.f95 ComputeD0.f95 ComputeDG82.f95 ComputeDm.f95 DHaj.f95 \
	djpi2.f95 FFTW3.f95 HilmGLQ.f95 MakeGrid2D.f95 \
	GLQGridCoord.f95 MakeGridGLQ.f95 MakeGridPoint.f95 \
	MakeMagGrid2D.f95 PlanetsConstants.f95 PlBar.f95 PlBar_d1.f95 \
	PLegendre.f95 PLegendreA.f95 PLegendreA_d1.f95 PlmBar.f95 \
	PlmBar_d1.f95 PlmIndex.f95 PlmSchmidt.f95 PlmSchmidt_d1.f95 \
	PlSchmidt.f95 PlSchmidt_d1.f95 PreCompute.f95 PreGLQ.f95 Random.f95 \
	ReadWisdom.f95 SHAdmitCorr.f95 SHBias.f95 \
	SHBiasK.f95 SHConvertCoef.f95 \
	SHExpandDH.f95 SHExpandGLQ.f95 \
	SHFindLWin.f95 SHLocalizedAdmitCorr.f95 SHMultiply.f95 \
	SHMultiTaperCSE.f95 SHMultiTaperSE.f95 SHPowerSpectra.f95 SHRead.f95 \
	SHReadJPL.f95 SHReturnTapersM.f95 SHRotateCoef.f95 \
	SHRotateRealCoef.f95 SHTOOLS.f95 SphericalCapCoef.f95 \
	Wigner3j.f95 wl.f95 SHReadGRACE.f95 SHReadCHAMP.f95 MakeGeoidGrid.f95 \
	MakeGridDH.f95 MakeCircleCoord.f95 SHSjkPG0.f95  \
	SHReturnTapers.f95 SHSjkPG.f95 PlON.f95 PlON_d1.f95 \
	PlmON.f95 PlmON_d1.f95 MakeGravGrid2D.f95 SHConfidence.f95 SHMagPowerSpectra.f95 \
	SHExpandDHC.f95 MakeGridDHC.f95 SHExpandGLQC.f95 MakeGridGLQC.f95 \
	SHPowerSpectraC.f95 SHBiasAdmitCorr.f95 SHCimToVector.f95 YilmIndex.f95 \
	ComputeDMap.f95 SHReturnTapersMap.f95 Curve2Mask.f95

OBJS0 =	CilmPlus.o ComputeD0.o ComputeDG82.o ComputeDm.o DHaj.o djpi2.o \
	FFTW3.o HilmGLQ.o MakeGrid2D.o GLQGridCoord.o MakeGridGLQ.o \
	MakeGridPoint.o MakeMagGrid2D.o PlanetsConstants.o PlBar.o PlBar_d1.o \
	PLegendre.o PLegendreA.o PLegendreA_d1.o PlmBar.o PlmBar_d1.o \
	PlmIndex.o PlmSchmidt.o PlmSchmidt_d1.o PlSchmidt.o PlSchmidt_d1.o \
	PreCompute.o PreGLQ.o Random.o ReadWisdom.o SHAdmitCorr.o \
	SHBias.o SHBiasK.o  \
	SHConvertCoef.o SHExpandDH.o SHExpandGLQ.o \
	SHFindLWin.o SHLocalizedAdmitCorr.o  \
	SHMultiply.o SHMultiTaperCSE.o SHMultiTaperSE.o SHPowerSpectra.o \
	SHRead.o SHReadJPL.o SHReturnTapersM.o SHRotateCoef.o \
	SHRotateRealCoef.o SHTOOLS.o SphericalCapCoef.o Wigner3j.o \
	wl.o SHReadGRACE.o SHReadCHAMP.o MakeGeoidGrid.o MakeGridDH.o \
	MakeCircleCoord.o SHSjkPG0.o SHReturnTapers.o SHSjkPG.o \
	PlON.o PlON_d1.o PlmON.o PlmON_d1.o MakeGravGrid2D.o \
	SHConfidence.o SHMagPowerSpectra.o SHExpandDHC.o MakeGridDHC.o \
	SHExpandGLQC.o MakeGridGLQC.o SHBiasAdmitCorr.o SHCilmToVector.o \
	YilmIndex.o ComputeDMap.o SHReturnTapersMap.o Curve2Mask.o
	
SRCS1 = EigValSym.f95 EigValVecSym.f95 EigValVecSymTri.f95 SHExpandLSQ.f95 SHMTDebias.f95  \
	SHMTVarOpt.f95 SHMTVarOpt0.f95

SRCS2 = EigValSym2.f95 EigValVecSym2.f95 EigValVecSymTri2.f95 SHExpandLSQ2.f95 SHMTDebias2.f950 \
	SHMTVarOpt2.f95 SHMTVarOpt02.f95

OBJS1 = EigValSym.o EigValVecSym.o EigValVecSymTri.o SHExpandLSQ.o SHMTDebias.o \
	SHMTVarOpt0.o SHMTVarOpt.o

OBJS2 = EigValSym2.o EigValVecSym2.o EigValVecSymTri2.o SHExpandLSQ2.o SHMTDebias2.o \
	SHMTVarOpt02.o SHMTVarOpt2.o


ifeq ($(MAKECMDGOALS),all)
	SRCS = $(SRCS0) $(SRCS1)
	OBJS = $(OBJS0) $(OBJS1)	
endif

ifeq ($(MAKECMDGOALS),all2)
	SRCS = $(SRCS0) $(SRCS2)
	OBJS = $(OBJS0) $(OBJS2)
endif


all: $(PROG) 

all2: $(PROG)

$(PROG): $(OBJS) 
	@echo
	@echo Compilation of source files successful.
	@echo
	@rm -f $(PROG)
#	@$(LIBTOOL) $(LIBTOOLFLAGS) -o $(PROG) $(OBJS)
#	If you prefer to use libtool, uncomment the above line, and comment the two lines below (AR and RLIB)
	-$(AR) $(ARFLAGS) $(PROG) $(OBJS) 
	-$(RLIB) $(RLIBFLAGS) $(PROG) 
	@echo
	@echo Creation of static library successful.
	@rm -f $(OBJS)
	@mv -f *.mod ../modules/
	-@cp -f *.mod ../modules/
	@echo Library and module files moved to ../lib and ../modules
	@echo
	@echo Archive successfully created


.PHONY: clean all all2

clean:
	rm -f $(PROG)
	rm -f $(OBJS0) 
	rm -f $(OBJS1) 
	rm -f $(OBJS2) 
	rm -f *.mod 
	rm -f ../modules/*.mod
	rm -f pyshtools.so pyshtools.o


.SUFFIXES: $(SUFFIXES) .f95 


.f95.o:
	$(F95) $(F95FLAGS) -c $<
	


CilmPlus.o: SHTOOLS.o
ComputeD0.o: SHTOOLS.o
ComputeDm.o: SHTOOLS.o
HilmGLQ.o: SHTOOLS.o
MakeGrid2D.o: SHTOOLS.o
GLQGridCoord.o: SHTOOLS.o
MakeGridGLQ.o: FFTW3.o SHTOOLS.o
MakeGridDH.o: FFTW3.o SHTOOLS.o
MakeGridPoint.o: SHTOOLS.o
MakeMagGrid2D.o: SHTOOLS.o
PreCompute.o: SHTOOLS.o
SHAdmitCorr.o: SHTOOLS.o
SHBias.o: SHTOOLS.o
SHBiasK.o: SHTOOLS.o
SHDeltaX.o: SHTOOLS.o
SHExpandDH.o: FFTW3.o SHTOOLS.o
SHExpandGLQ.o: FFTW3.o SHTOOLS.o
SHExpandLSQ.o: SHTOOLS.o
SHExpandLSQ2.o: SHTOOLS.o
SHFindLWin.o: SHTOOLS.o
SHLocalizedAdmitCorr.o: SHTOOLS.o
SHMultiply.o: SHTOOLS.o
SHMultiTaperCSE.o: SHTOOLS.o
SHMultiTaperSE.o: SHTOOLS.o
SHReturnTapersM.o: SHTOOLS.o
SHReturnTapers.o : SHTOOLS.o
SHRotateRealCoef.o: SHTOOLS.o
SHSjkPG0.o: SHTOOLS.o
SHSjkPG.o: SHTOOLS.o
SphericalCapCoef.o: SHTOOLS.o
MakeGeoidGrid.o: SHTOOLS.o
SHMTVarOpt0.o : SHTOOLS.o
SHMTVarOpt.o : SHTOOLS.o
PLegendreA.o  : SHTOOLS.o
PLegendreA_d1.o : SHTOOLS.o
PlmBar.o : SHTOOLS.o
PlmBar_d1.o : SHTOOLS.o
PlmSchmidt.o : SHTOOLS.o
PlmSchmidt_d1.o : SHTOOLS.o
PlmON.0 : SHTOOLS.o
PlmON_d1.o : SHTOOLS.o
SHMTDebias.o : SHTOOLS.o
SHMTDebias2.o : SHTOOLS.o
MakeGravGrid2d.o : SHTOOLS.o
SHBiasAdmitCorr.o : SHTOOLS.o
ComputeDMap.o : SHTOOLS.o
SHReturnTapersMap.o : SHTOOLS.o
