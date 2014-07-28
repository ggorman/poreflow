Poreflow
========

Modelling flow at the rock pore scale

Compilation
---------
```
cd src/
make
```
The binaries are stored in bin/

Hourglass example
-----------------
In this simple example we:
* Create a simple image of a pore channel with a hourglass profile.
* Generate a tetrahedral mesh using the Tarantula meshing code.
* Use dolfin to model stokes flow through the channel and calculate the flux and permeability.

```
# Create an image 47^3 voxels in size, norrowest width of hourglass. Write file in vox image format.
./bin/create_hourglass -s 47 -t 8 -c vox -o hourglass.vox

# Mesh image using tarantula's voxmesher
voxmesher hourglass.conf

# Convert Tarantula's .spm file into GMSH format.
./bin/tarantula2gmsh -v hourglass.spm

# Convert GMSH file into Dolfin's xml format.
dolfin-convert hourglass.msh hourglass.xml

# Run Stokes FEM model using a direct solver and the Taylor-Hood element pair.
python python/stokes-dolfin.py -D -e 0 hourglass.xml 
```
You will see that the file velocity000000.vtu has been created. Use ParaView to visualise the velocity profile.

Small Berea Sandstone example
-----------------------------
In this simple example we:
* Extract 64x64x64 section out of the Berea sandstone sample from http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling/micro-ct%20images%20and%20networks/berea%20sandstone
* Generate a tetrahedral mesh using the Tarantula meshing code.
* Use dolfin to model stokes flow through the channel and calculate the flux and permeability.

```
# Extract 64^3 voxels from micro-CT image. Write out file in vox image format.
cd data/
unzip 33505696.ZIP 
mv Image Berea
./bin/convert_microct -c vox -s 64 Berea

# Mesh image using tarantula's voxmesher
voxmesher Berea.conf

# Convert Tarantula's .spm file into GMSH format.
../bin/tarantula2gmsh -v Berea.spm

# Convert GMSH file into Dolfin's xml format.
dolfin-convert Berea.msh Berea.xml

# Run Stokes FEM model using a direct solver and the Taylor-Hood element pair.
python python/stokes-dolfin.py -D -e 4 Berea.xml
```
You will see that the file velocity000000.vtu has been created. Use ParaView to visualise the velocity profile.

