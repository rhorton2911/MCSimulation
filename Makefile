#Object files to produce and link together
e = intp.o numbers.o input_data.o states.o elasticdcs.o totalcs.o tcstostate.o sdcs.o PsFormation.o montecarlo.o AnalyticScattering.o ParticleDynamics.o simulation.o main.o

###-------------------------------ifort compiler setup------------------------------###
###DEBUGGING###
#a = -O0 -C -g  -r8 -shared-intel -qopenmp -traceback -heap-array -dynamic 
#-mcmodel=medium
#DEBUGGING flags for use with ARM profiling tools(requires compiler optimisation to be turned off and static linking for tools to work)#a = -O0 -C -g -r8 -qopenmp -traceback -heap-array -dynamic 
#a = -O0 -C -g -r8 -qopenmp -traceback -heap-array -dynamic 
###FURTHER DEBUGGING###
#a := $(a) -debug extended -check all -warn all -ftrapuv -fp-stack-check
#Enable all optimisations
#a = -O3 -heap-array -dynamic
#CFLAGS := -O3 -xHost -ip -fpp -qopenmp -no-wrap-margin -fp-model precise -prec-div \
            -prec-sqrt -qoverride-limits -traceback -reentrancy threaded \
	    -mcmodel=medium -shared-intel -r8 -align dcommons \
	    -h PIC -diag-disable 8291 -parallel  #-mkl
#a := $(a) $(CFLAGS)
#l:= -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

#Recommended link line from intel link line advisor
#l:= -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

#Modify link line when using arm ddt debugging tools, static link to ddt libraries
#l:= -Wl, --wrap=dlopen,--wrap=dlclose,--allow-multiple-definition,--undefined=malloc/$(LD_LIBRARY_PATH)/64/libdmallocth.a $(l)
#!!! current use
# a = -O3  -r8 

# O0 : disables all optimisations.
# C : checks for certain conditions at run-time, including temporary arguments, assume directives, indexing bounds, format types and lengths, pointers, stacks, and variable initialisation. For more info search -check in man ifort.
# g : Generate full debugging information in the object file.
# openmp : use OpenMP.
# traceback : Provide extra information when a severe error occurs at run time.
# r8 : Default size for real (complex) variables are *8 (*16) bytes long.
# fpp : Runs the Fortran preprocessor on source files before compilation.



##------------------------------------------- gfortran compiler setup -------------------------------##
###DEBUGGING FLAGS###
a = -Wall -Wextra -pedantic -Warray-bounds -fbacktrace -fopenmp -frecursive -g -Og \
    -fimplicit-none -fcheck=all -ffree-line-length-none -fdefault-integer-8 -m64 -Wno-tabs

#fclink = -fdefault-integer-8 -m64 

###Production Flags - Enable Optimisations###
#a = -ffree-line-length-none -fopenmp -fdefault-integer-8 -m64
##---------------------------------------------------------------------------------------------------##
 
#ftn is a fortran compiler wrapper provided by Pawsey. Compiler used is selected based on the current programming environment.
#To use intel, load PrgEnv-Intel, to use gfortran load PrgEnv-gnu 
#Wrapper handles linking to LAPACK libraries when using gnu environment, no extra link line commands required.
FC = ftn  

export a
export PRECISION

main: $e 
	$(FC) $a $e  -o main

clean:
	rm -f *.mod
	rm -f *.o

#Must be compiled first, other modules require data types in this one
numbers.o : numbers.f90
	$(FC) $a -c numbers.f90

input_data.o : input_data.f90 numbers.f90
	$(FC) $a -c input_data.f90

main.o :  main.f90 numbers.f90 input_data.f90 sdcs.f90 states.f90 totalcs.f90 tcstostate.f90
	$(FC) $a -c main.f90

states.o : states.f90 numbers.f90 input_data.f90 
	$(FC) $a -c states.f90

totalcs.o : totalcs.f90 numbers.f90 input_data.f90 montecarlo.f90
	$(FC) $a -c totalcs.f90

tcstostate.o : tcstostate.f90 numbers.f90 input_data.f90 
	$(FC) $a -c tcstostate.f90

sdcs.o : sdcs.f90 numbers.f90 input_data.f90 states.f90
	$(FC) $a -c sdcs.f90

intp.o : intp.f
	$(FC) $a -w -std=legacy -c intp.f

montecarlo.o : montecarlo.f90 numbers.f90 input_data.f90 totalcs.f90 states.f90
	$(FC) $a -c montecarlo.f90

AnalyticScattering.o : AnalyticScattering.f90 numbers.f90 montecarlo.f90 input_data.f90 totalcs.f90 states.f90 simulation.f90
	$(FC) $a -c AnalyticScattering.f90

ParticleDynamics.o : ParticleDynamics.f90 numbers.f90  montecarlo.f90 input_data.f90 totalcs.f90 states.f90 simulation.f90 AnalyticScattering.f90
	$(FC) $a -c ParticleDynamics.f90

simulation.o : simulation.f90 numbers.f90  montecarlo.f90 input_data.f90 totalcs.f90 states.f90 ParticleDynamics.f90 AnalyticScattering.f90
	$(FC) $a -c simulation.f90

PsFormation.o : PsFormation.f90 numbers.f90 simulation.f90 montecarlo.f90
	$(FC) $a -c PsFormation.f90

elasticdcs.o : elasticdcs.f90 numbers.f90 input_data.f90
	$(FC) $a -c elasticdcs.f90


