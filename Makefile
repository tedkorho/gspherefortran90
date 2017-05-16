LOWLEVELMODS = 	lib/specfunc.f90 lib/randev.f90 lib/voper.f90
MODS = 		lib/corrfunc.f90 lib/discrets.f90 lib/gsaxg.f90 lib/gsg.f90
MAINPROG =	lib/gsphere.f90
OBJS =		gsphere.o discrets.o corrfunc.o gsg.o gsaxg.o randev.o specfunc.o voper.o
MODFILES =	discrets.mod corrfunc.mod gsaxg.mod gsg.mod randev.mod specfunc.mod voper.mod

all: main

clean:
	rm $(OBJS) $(MODFILES)

main: lib/gsphere.f90 lib/corrfunc.f90 lib/discrets.f90 lib/gsg.f90 lib/gsaxg.f90 lib/randev.f90 lib/specfunc.f90 lib/voper.f90
	gfortran -c $(LOWLEVELMODS)
	gfortran -c $(MODS)
	gfortran -c $(MAINPROG)
	gfortran -o gsphere $(OBJS)
	rm -f *.o
	rm -f *.mod
