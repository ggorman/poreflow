# Solves the Stokes equation. This is adapted from the fenics example:
# http://fenicsproject.org/documentation/dolfin/dev/python/demo/pde/stokes-iterative/python/documentation.html
#
# Note that the sign for the pressure has been flipped for symmetry."""

from dolfin import *
from numpy import array
from petsc4py import PETSc

import getopt

import sys

parameters["form_compiler"]["quadrature_degree"] = 6
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["cpp_optimize_flags"] = "-O3 -ffast-math -march=native"

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
      -w L  Zidth of domain.
"""

# Test for PETSc
if not has_linear_algebra_backend("PETSc"):
    info("DOLFIN has not been configured with PETSc. Exiting.")
    exit()

optlist, args = getopt.getopt(sys.argv[1:], 'dhve:')

verbose = False
element_pair = 4
sample_width = None
direct_solver = False
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
  elif opy[0] == '-w':
    sample_width = float(opt[1])  

filename = args[-1]

try:
  axis = int(filename[:-4].split('_')[-4])
  bc_in = int(filename[:-4].split('_')[-2])
  bc_out = int(filename[:-4].split('_')[-1])
except:
  axis = 0
  bc_in = 1
  bc_out = 2

lu_args = [sys.argv[0]] + """
                             --petsc.se_snes_monitor
                             --petsc.se_snes_converged_reason
                             --petsc.se_snes_stol 0.0
                             --petsc.se_snes_atol 1.0e-9
                             --petsc.se_snes_rtol 1.0e-9
                             --petsc.se_snes_max_it 200
                             --petsc.se_ksp_type preonly
                             --petsc.se_pc_type lu
                             --petsc.se_pc_factor_mat_solver_package mumps
                             """.split()

fieldsplit_args = [sys.argv[0]] + """
                             --petsc.se_snes_monitor
                             --petsc.se_snes_max_it 30
                             --petsc.se_snes_converged_reason
                             --petsc.se_snes_rtol 1.0e-7
                             --petsc.se_snes_atol 1.0e-7

                             --petsc.se_ksp_converged_reason
                             --petsc.se_ksp_type gcr
                             --petsc.se_ksp_monitor_true_residual
                             --petsc.se_ksp_rtol 1.0e-7
                             --petsc.se_ksp_atol 1.0e-7

                             --petsc.se_pc_type fieldsplit
                             --petsc.se_pc_fieldsplit_type schur
                             --petsc.se_pc_fieldsplit_schur_factorization_type upper

                             --petsc.se_fieldsplit_0_ksp_type preonly
                             --petsc.se_fieldsplit_0_ksp_max_it 1
                             --petsc.se_fieldsplit_0_pc_type ml

                             --petsc.se_fieldsplit_1_ksp_type preonly
                             --petsc.se_fieldsplit_1_ksp_max_it 2
                             --petsc.se_fieldsplit_1_pc_type jacobi
                        """.split()

parameters["std_out_all_processes"] = False
args = fieldsplit_args
parameters.parse(args)

mesh = Mesh(filename)

n = FacetNormal(mesh)

# Define function spaces 
Z = None
if element_pair == 0:
  # - Taylor-Hood
  p = 2
  V = VectorFunctionSpace(mesh, "Lagrange", p)
  Q = FunctionSpace(mesh, "Lagrange", p-1)
  Z = V * Q
elif element_pair == 1:
  V = VectorFunctionSpace(mesh, "Lagrange", 3)
  Q = FunctionSpace(mesh, "Lagrange", 1)
  Z = V * Q
elif element_pair == 2:
  # - Crouzeix-Raviart
  V = VectorFunctionSpace(mesh, "Crouzeix-Raviart", 1)
  Q = FunctionSpace(mesh, "Discontinuous Lagrange", 0)
  Z = V * Q
elif element_pair == 3:
  # - CD
  V = VectorFunctionSpace(mesh, "Lagrange", 3)
  Q = FunctionSpace(mesh, "Discontinuous Lagrange", 1)
  Z = V * Q
elif element_pair == 4:
  # - MINI
  p = 1
  P1 = VectorFunctionSpace(mesh, "Lagrange", p)
  B = VectorFunctionSpace(mesh, "Bubble", p+3)
  Q = FunctionSpace(mesh, "Lagrange", p)
  Z = (P1+B)*Q

# Boundary
boundaries = MeshFunction('size_t', mesh, filename[0:-4] + "_facet_region.xml")

ds = Measure('ds')[boundaries]

# A few physical paramerters
dP = 1.0

# Boundary conditions
noslip = Constant((0.0, 0.0, 0.0))
p_in = Constant(dP)
p_out = Constant(0.0)
bcs = []
for i in range(1, 8):
    if i != bc_in and i != bc_out:
        bcs.append(DirichletBC(Z.sub(0), noslip, boundaries, i))

# Define variational problem
z      = Function(Z)
(u, p) = split(z)
(v, q) = TestFunctions(Z)

F = (inner(grad(u), grad(v)) - div(v)*p + q*div(u))*dx + p_in*dot(v, n)*ds(bc_in)

# Class for interacting with Newton solver.
class GeneralProblem(NonlinearProblem):
  def __init__(self, F, z, bcs):
    NonlinearProblem.__init__(self)
    self.fform = F
    self.z = z
    self.bcs = bcs
    self.jacobian = derivative(F, z)

  def F(self, b, x):
    assemble(self.fform, tensor=b)
    [bc.apply(b, x) for bc in self.bcs]

  def J(self, A, x):
    assemble(self.jacobian, tensor=A)
    [bc.apply(A) for bc in self.bcs]

problem = GeneralProblem(F, z, bcs=bcs)
solver  = PETScSNESSolver()

solver.parameters["relative_tolerance"] = 1e-6
solver.parameters["report"] = False
solver.parameters["options_prefix"] = "se"

solver.init(problem, z.vector())
snes = solver.snes()
snes.setFromOptions()

# Configure the FIELDSPLIT stuff.
if args != lu_args:
  pc = snes.ksp.pc

  fields = []
  for i in range(2):
    subspace = SubSpace(Z, i)
    subdofs  = subspace.dofmap().dofs()
    IS = PETSc.IS()
    IS.createGeneral(subdofs)
    name = str(i)
    fields.append((name, IS))

  pc.setFieldSplitIS(*fields)

  def extract_sub_matrix(A, subspace_in, subspace_out):
    Amat  = as_backend_type(A).mat()

    subis_in   = PETSc.IS()
    subdofs_in = Z.sub(subspace_in).dofmap().dofs()
    subis_in.createGeneral(subdofs_in)

    subis_out   = PETSc.IS()
    subdofs_out = Z.sub(subspace_out).dofmap().dofs()
    subis_out.createGeneral(subdofs_out)

    submat  = Amat.getSubMatrix(subis_out, subis_in)
 
    return submat

  def extract_sub_vector(V, subspace):
    Vvec  = as_backend_type(V.vector()).vec()

    subis   = PETSc.IS()
    subdofs = Z.sub(subspace).dofmap().dofs()
    subis.createGeneral(subdofs)

    subvec  = Vvec.getSubVector(subis)
    dupe    = subvec.duplicate()
    Vvec.restoreSubVector(subis, subvec)
 
    return dupe

  (u, p) = TrialFunctions(Z)
  (v, q) = TestFunctions(Z)

  schur_D = assemble(inner(p, q)*dx)
  #[bc.apply(schur_D) for bc in bcs]
  schur = extract_sub_matrix(schur_D, 1, 1)

  def monitor(snes, its, norm):
    pc = snes.ksp.pc
    pc.setFieldSplitSchurPreType(PETSc.PC.SchurPreType.USER, schur)

  snes.setMonitor(monitor)

snes.solve(None, as_backend_type(z.vector()).vec())

File("solution_stokes.xml.gz") << z

(u, p) = z.split()

# Darcy's law
# U = -k/mu dP/L
# Darcy velocity U = flux/area = flux/L^2
# Therefore, permability, k = UL/dP = (flux mu)/(L dP) (here mu==1)

flux = [assemble(dot(u, n)*ds(i)) for i in range(1, 8)]
mean_flux = (-flux[bc_in-1]+flux[bc_out-1])*0.5

bbox = array([MPI.min(mesh.mpi_comm(), mesh.coordinates()[:,0].min()), MPI.max(mesh.mpi_comm(), mesh.coordinates()[:,0].max()),
              MPI.min(mesh.mpi_comm(), mesh.coordinates()[:,1].min()), MPI.max(mesh.mpi_comm(), mesh.coordinates()[:,1].max()),
              MPI.min(mesh.mpi_comm(), mesh.coordinates()[:,2].min()), MPI.max(mesh.mpi_comm(), mesh.coordinates()[:,2].max())])

L = (bbox[1]-bbox[0] + bbox[3]-bbox[2] + bbox[5]-bbox[4])/3.0

if sample_width:
    mean_flux*=(sample_width/L)**3
    L = sample_width

permability = (mean_flux)/(dP*L)

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

