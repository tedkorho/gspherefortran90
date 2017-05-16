main: lib/gsphere.f90 lib/corrfunc.f90 lib/discrets.f90 lib/gsg.f90 lib/gsaxg.f90 lib/randev.f90 lib/specfunc.f90 lib/voper.f90
	gfortran -c lib/corrfunc.f90 lib/discrets.f90 lib/gsg.f90 lib/gsaxg.f90 lib/randev.f90 lib/specfunc.f90 lib/voper.f90
	gfortran -c gsphere.f90
	gfortran -o gsphere gsphere.o discrets.o corrfunc.o gsg.o gsaxg.o randev.o specfunc.o voper.o
