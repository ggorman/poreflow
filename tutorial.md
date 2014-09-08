Tutorial for the impatient
==========================

Download and compile
--------------------
If you are working in AMCG/ESE then do these initial steps on valhalla.ese.ic.ac.uk as that is where the Tarantula mesh generator licensed.

```bash
mkdir ~/projects/
cd projects
git clone https://github.com/ggorman/poreflow.git
cd poreflow/src
make
echo "export PATH=$HOME/projects/poreflow/bin:$PATH" >> .bashrc
export PATH=$HOME/projects/poreflow/bin:$PATH
```

Prepare data sample
-------------------
All the data is available [here][http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling/micro-ct%20images%20and%20networks]. For this example click on [Berea Sandstone][http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling/micro-ct%20images%20and%20networks/berea%20sandstone]. You can see a full description of the data sample. Download the image by either clicking on [Download image (binary and ASCII)][http://www3.imperial.ac.uk/pls/portallive/docs/1/33505696.ZIP] or it may be easier to use wget:
```bash
wget http://www3.imperial.ac.uk/pls/portallive/docs/1/33505696.ZIP
unzip 33505696.ZIP
mv Image ~/projects/Berea
```

Next you have to convert the raw format to something that the mesh generator can read (there is going to be a lot of converting so get used to it).
```bash
convert_microct -v -c vox -s 64 Berea
```

If you just type the command by itself it will give you an example of its usage. The *-v* option enables verbose mode. As well as outputting messages to the stdout it may also output VTK files for diagnostic/visualisation purposes. *-c vox* specifies the output format as VOX and the optional argument *-s <size>* specifies that only a sample of the data is selected.

These files can be huge so always compress them after generation or you will quickly run out of space and patience.

```bash
gzip Berea.vox
```

Mesh generation
---------------
Grab a template configuration file for Tarantula:

```bash
cp $HOME/projects/poreflow/data/Berea.conf .
```

Although obviously not required here, rename this to the rename this file to correspond to the dataset you are working on, and update the first 2 lines of the conf file so that the vox and mesh names are correct. I strongly recommend you use a different directory for each data sample or you will continuously clobber old results/data. Next, run the mesh generator:

```bash
/usr/local/tarantula/bin/linux/voxmesher Berea.conf
```

This can sometimes fail. When this happens try modifying the sample size with the -s option above. Powers of 2 are good. If this runs successfully you should see a Tarantula .spm file:

```bash
ls -ltr | tail -1 
-rw-r--r--  1 ggorman  esestaff   65685433 Sep  8 11:23 Berea.spm
```

This us a rather unique format that not many codes will recognise so we will convert this to GMSH format. For this we will use tarantula2gmsh. This will also just extract the void region.

```bash
tarantula2gmsh -t -v Berea.spm
```

As well as converting the file format, this command does the following:

* Images will typically have two materials (rock and void indicated by 1 and 0), you can toggle which material mesh it extracts using the *-t* flag.
* Extracts only the active region - it throws away any connected region that is not connected to both sides of the domain along the X-axis.
* Applies boundary labels: -x, +x, -y, +y, -z, +z, grain boundaries labelled as 1, 2, 3, 4, 5, 6, 7 respectively.
* Add the *-v* option if you want verbose messaging and VTK files to admire your beautiful mesh!

Use paraview to take a look at the data. Does it look ok? Is it "fit for purpose"?

Running a simulation
--------------------
We are going to use caloris.ese.ic.ac.uk because this has the master version of FEniCS and PETSc installed (complements of Patrick Farrell) which is required for the split field preconditioners used.

First set environment variables so you are using the right versions of everything:

```bash
source /data/pfarrell/src/local/install_fenics_opt.sh
```

Convert GMSH to Dolfin XML format so that dolfin can read it. Ideally we should convert to XDMF but I am having some trouble getting that working:

```bash
dolfin-convert Berea.msh Berea.xml
python ~/projects/poreflow/python/convert_dolfin_new_xml.py Berea.xml
```

Finally we are ready to run the Stokes flow solver in parallel using MPI:

```bash
mpiexec -n 24 python ~/projects/poreflow/python/stokes-dolfin-snes.py Berea.xml

 # the tail of the output will look something like this:
 176 KSP preconditioned resid norm 1.074899229318e-05 true resid norm 2.028839668985e-06 ||r(i)||/||b|| 6.432845983965e-08
 177 KSP preconditioned resid norm 1.078483993909e-05 true resid norm 2.035275971623e-06 ||r(i)||/||b|| 6.453253581576e-08
 178 KSP preconditioned resid norm 5.474046571645e-06 true resid norm 1.032561103527e-06 ||r(i)||/||b|| 3.273943549884e-08
 Linear solve converged due to CONVERGED_RTOL iterations 178
 1 SNES Function norm 1.032561103524e-06 
 Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 1
 [-0.7067218987108811, 0.7067164213390983, 0.0, 0.0, 0.0, 0.0, 0.0]
 #################################
 ## In-flux = -0.706722
 ## Out-flux = 0.706716
 ## Permability = 0.0112153
 #################################
```

Note the permability above. Does it look the value you were expecting?

Finally have a look at the solution using paraview:

```bash
paraview --data=velocity.pvd 
paraview --data=pressure.pvd 
```
