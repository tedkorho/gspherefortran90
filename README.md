GSPHERE-FORTRAN 90
Updated to modern Fortran spec from original code
by Karri Muinonen.
Teo Korhonen

Generates a Gaussian sphere with randomly distributed surface height.
Uses triangle discretization to create a 3D model. 

COMPILATION AND EXECUTION:

Requires gfortran.
Compile with the command

> make

in the program's root folder.

The executable is ./gsphere. Requires an input file
argument formatted like the example files in input folder.

OUTPUT AND VISUALIZATION:

Outputs in 3 formats: Matlab X,Y, and Z-coordinates, .idf file, 
and .vtk file.

Matlab X,Y,Z:

Use any suitable f

> float = '%f'
> xfile = output/mloutx; yfile = output/mlouty; zfile = output/mloutz
> 
> x = fscanf(xfile,float); y = fscanf(yfile,float); z = fscanf(zfile,float)
> 
> scatter3(x,y,z)

.idf: Use any compatible visualization tool.

.vtk: ParaView, visit or other compatible visualization tools can 
open it right away.

OTHER:

test.py script requires Valgrind.
