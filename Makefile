PARAM =		src/parameters.f90
LOWLEVELMODS = 	src/specfunc.f90 src/randev.f90 src/voper.f90
MODS = 		src/corrfunc.f90 src/discrets.f90 src/gsaxg.f90 src/gsg.f90
MAINPROG =	src/gsphere.f90
OBJS =		gsphere.o discrets.o corrfunc.o gsg.o gsaxg.o randev.o specfunc.o voper.o
MODFILES =	discrets.mod corrfunc.mod gsaxg.mod gsg.mod randev.mod specfunc.mod voper.mod

all: main

main: $(MAINPROG) $(MODS) $(LOWLEVELMODS) $(PARAM)
	gfortran -c $(PARAM)
	gfortran -c $(LOWLEVELMODS) -Wall
	gfortran -c $(MODS) -Wall
	gfortran -c $(MAINPROG) -Wall
	gfortran -o gsphere $(OBJS) -Wall
	rm -f *.o
	rm -f *.mod
