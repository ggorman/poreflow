# Solves the Stokes equations using an iterative linear solver.
# Taken from:
# http://fenicsproject.org/documentation/dolfin/dev/python/demo/pde/stokes-iterative/python/documentation.html
#
# Note that the sign for the pressure has been flipped for symmetry."""

from dolfin import *
from numpy import array

import sys

set_log_level(DEBUG)

# Test for PETSc or Epetra
if not has_linear_algebra_backend("PETSc") and not has_linear_algebra_backend("Epetra"):
    info("DOLFIN has not been configured with Trilinos or PETSc. Exiting.")
    exit()

if not has_krylov_solver_preconditioner("amg"):
    info("Sorry, this demo is only available when DOLFIN is compiled with AMG "
	 "preconditioner, Hypre or ML.");
    exit()

filename = sys.argv[1]
try:
    sample_width = float(sys.argv[2])
except:
    sample_width = None

try:
	axis = int(filename[:-4].split('_')[-4])
	bc_in = int(filename[:-4].split('_')[-2])
	bc_out = int(filename[:-4].split('_')[-1])
except:
	axis = 0
	bc_in = 1
	bc_out = 2

mesh = Mesh(filename)

n = FacetNormal(mesh)


# Define function spaces 
# - Taylor-Hood
#p = 2
#V = VectorFunctionSpace(mesh, "Lagrange", p)
#Q = FunctionSpace(mesh, "Lagrange", p-1)
#W = V * Q

# V = VectorFunctionSpace(mesh, "Lagrange", 3)
# Q = FunctionSpace(mesh, "Lagrange", 1)
# W = V * Q

# - Crouzeix-Raviart
# V = VectorFunctionSpace(mesh, "Crouzeix-Raviart", 1)
# Q = FunctionSpace(mesh, "Discontinuous Lagrange", 0)
# W = V * Q

# - CD
# V = VectorFunctionSpace(mesh, "Lagrange", 3)
# Q = FunctionSpace(mesh, "Discontinuous Lagrange", 1)
# W = V * Q

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
noslip = Constant((0.0, 0.0, 0.0))
p_in = Constant(-1.0)
p_out = Constant(0.0)
bcs = []
for i in range(1, 8):
    if i == bc_in: # Inflow
        bcs.append(DirichletBC(W.sub(1), p_in, boundaries, bc_in))   
    elif i == bc_out: # Outflow
        bcs.append(DirichletBC(W.sub(1), p_out, boundaries, bc_out))
    else:
        bcs.append(DirichletBC(W.sub(0), noslip, boundaries, i))

# Define variational problem
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)
f = Constant((0.0, 0.0, 0.0))
# a = inner(grad(u), grad(v))*dx + div(v)*p*dx + q*div(u)*dx - (p*dot(v, n)*ds(bc_in) + p*dot(v, n)*ds(bc_out))
a = inner(grad(u), grad(v))*dx + div(v)*p*dx + q*div(u)*dx
L = inner(f, v)*dx

direct_solver = True

# Form for use in constructing preconditioner matrix
b = inner(grad(u), grad(v))*dx + p*q*dx

# Assemble system
A, bb = assemble_system(a, L, bcs)

if direct_solver:
    solver = LUSolver("mumps")
    # solver = LUSolver("umfpack")
    solver.set_operator(A)
else:
    # Assemble preconditioner system
    P, btmp = assemble_system(b, L, bcs)

    # Create Krylov solver and AMG preconditioner
    # solver = KrylovSolver("tfqmr", "amg")
    solver = KrylovSolver("gmres", "amg")
    solver.parameters['absolute_tolerance'] = 1E-6
    solver.parameters['relative_tolerance'] = 1E-6
    solver.parameters['maximum_iterations'] = 100
    solver.parameters['monitor_convergence'] = True
    solver.parameters['report'] = True
    solver.parameters['error_on_nonconvergence'] = False

    # Associate operator (A) and preconditioner matrix (P)
    solver.set_operators(A, P)

# Solve
U = Function(W)
solver.solve(U.vector(), bb)

# Get sub-functions
u, p = U.split()

# Darcy's law
# U = -k/mu dP/L
# Here the visocity, mu, and dP are both set to unity. The Darcy velocity U = flux/area = flux/L^2
# Therefore, permability, k = UL = flux/L

flux = [assemble(dot(u, n)*ds(i)) for i in range(1, 8)]
mean_flux = (-flux[bc_in-1]+flux[bc_out-1])*0.5

bbox = array([MPI.min(mesh.coordinates()[:,0].min()), MPI.max(mesh.coordinates()[:,0].max()),
              MPI.min(mesh.coordinates()[:,1].min()), MPI.max(mesh.coordinates()[:,1].max()),
              MPI.min(mesh.coordinates()[:,2].min()), MPI.max(mesh.coordinates()[:,2].max())])

L = (bbox[1]-bbox[0] + bbox[3]-bbox[2] + bbox[5]-bbox[4])/3.0

if sample_width:
    mean_flux*=(sample_width/L)**3
    L = sample_width

permability = mean_flux/L

if MPI.process_number() == 0:
    print flux
    print """#################################
## In-flux = %g
## Out-flux = %g
## Permability = %g
#################################"""%(flux[bc_in-1], flux[bc_out-1], permability)

# Because the sign of P was flipped
pflip = project(-p)

# Save solution in VTK format
ufile_pvd = File("velocity.pvd")
ufile_pvd << u
pfile_pvd = File("pressure.pvd")
pfile_pvd << pflip

# Plot solution
# plot(u)
# plot(pflip)
# interactive()
