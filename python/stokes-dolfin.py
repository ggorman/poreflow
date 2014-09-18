# Solves the Stokes equation. This is adapted from the fenics example:
# http://fenicsproject.org/documentation/dolfin/dev/python/demo/pde/stokes-iterative/python/documentation.html

from dolfin import *
from numpy import array

import getopt
import sys
import os.path

def usage():
  print sys.argv[0]+""" [options] dolfin_mesh.xml
    options:
      -h    Prints this help message.
      -d    Enable debugging mode.
      -v    Enable verbose mode.
      -e 0-4 where:
         0  P2 for velocity and P1 for pressure, i.e. Taylor-Hood.
         1  p3-p1.
         2  Crouzeix-Raviart.
         3  p3-DG P1, CD.
         4  MINI (default)
"""

# Test for PETSc
if not has_linear_algebra_backend("PETSc"):
    info("DOLFIN has not been configured with PETSc. Exiting.")
    exit()

# Set defaults for options.
verbose = False
element_pair = 4

# Parse commandline options.
optlist, args = getopt.getopt(sys.argv[1:], 'dhve:')
for opt in optlist:
  if opt[0] == '-d':
    set_log_level(DEBUG)
  elif opt[0] == '-h':
    usage()
    sys.exit(0)
  elif opt[0] == '-v':
    verbose = True
  elif opt[0] == '-e':
    element_pair = int(opt[1])

filename = args[-1]
if not os.path.isfile(filename):
    raise IOError("No such file: %s"%filename)
    usage()
    exit()

axis = 0
bc_in = 1
bc_out = 2

mesh = Mesh(filename)

n = FacetNormal(mesh)

# Define function spaces. 
W = None
if element_pair == 0:
  # - Taylor-Hood
  p = 2
  V = VectorFunctionSpace(mesh, "Lagrange", p)
  Q = FunctionSpace(mesh, "Lagrange", p-1)
  W = V * Q
elif element_pair == 1:
  V = VectorFunctionSpace(mesh, "Lagrange", 3)
  Q = FunctionSpace(mesh, "Lagrange", 1)
  W = V * Q
elif element_pair == 2:
  # - Crouzeix-Raviart
  V = VectorFunctionSpace(mesh, "Crouzeix-Raviart", 1)
  Q = FunctionSpace(mesh, "Discontinuous Lagrange", 0)
  W = V * Q
elif element_pair == 3:
  # - CD
  V = VectorFunctionSpace(mesh, "Lagrange", 3)
  Q = FunctionSpace(mesh, "Discontinuous Lagrange", 1)
  W = V * Q
elif element_pair == 4:
  # - MINI
  p = 1
  P1 = VectorFunctionSpace(mesh, "Lagrange", p)
  B = VectorFunctionSpace(mesh, "Bubble", p+3)
  Q = FunctionSpace(mesh, "Lagrange", p)
  W = (P1+B)*Q

# Boundary
boundaries = MeshFunction('size_t', mesh, filename[0:-4] + "_facet_region.xml")

ds = Measure('ds')[boundaries]

# Boundary conditions
dP = 1.0
noslip = Constant((0.0, 0.0, 0.0))
p_in = Constant(dP)
p_out = Constant(0.0)
bcs = []
for i in range(1, 8):
    if i != bc_in and i != bc_out:
        bcs.append(DirichletBC(W.sub(0), noslip, boundaries, i))

# Define variational problem
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

a = (inner(grad(u), grad(v)) - div(v)*p + q*div(u))*dx
A = assemble(a)

L = inner(Constant((0.0, 0.0, 0.0)), v)*dx - p_in*dot(v, n)*ds(bc_in)
b = assemble(L)

for bc in bcs:
    bc.apply(A)
    bc.apply(b)

solver = LUSolver("mumps")
solver.set_operator(A)

# Solve
U = Function(W)
solver.solve(U.vector(), b)

# Get sub-functions
u, p = U.split()

# Darcy's law
# U = -k/mu dP/L
# Here the visocity, mu, and dP are both set to unity. The Darcy velocity U = flux/area = flux/L^2
# Therefore, permability, k = UL = flux/L

flux = [assemble(dot(u, n)*ds(i)) for i in range(1, 8)]
mean_flux = (-flux[bc_in-1]+flux[bc_out-1])*0.5

bbox = array([MPI.min(mesh.mpi_comm(), mesh.coordinates()[:,0].min()), MPI.max(mesh.mpi_comm(), mesh.coordinates()[:,0].max()),
              MPI.min(mesh.mpi_comm(), mesh.coordinates()[:,1].min()), MPI.max(mesh.mpi_comm(), mesh.coordinates()[:,1].max()),
              MPI.min(mesh.mpi_comm(), mesh.coordinates()[:,2].min()), MPI.max(mesh.mpi_comm(), mesh.coordinates()[:,2].max())])

L = (bbox[1]-bbox[0] + bbox[3]-bbox[2] + bbox[5]-bbox[4])/3.0

permability = mean_flux/L

if MPI.rank(mesh.mpi_comm()) == 0:
    print flux
    print """#################################
## In-flux = %g
## Out-flux = %g
## Permability = %g
#################################"""%(flux[bc_in-1], flux[bc_out-1], permability)

# Save solution in VTK format
ufile_pvd = File("velocity.pvd")
ufile_pvd << u
pfile_pvd = File("pressure.pvd")
pfile_pvd << p

