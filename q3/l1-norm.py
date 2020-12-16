import numpy
from scipy import linalg
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 11
nx = 41  # number of spatial discrete points
L = 2.0  # length of the domain
dx = L / (nx - 1)  # spatial grid size
nt = 10  # number of time steps

# Define the grid point coordinates.
x = numpy.linspace(0.0, L, num=nx)
u0 = numpy.ones(nx) # base initial condition with u=1
# Get a list of indices where 0.5 >0 where u=2 as initial condition.
mask = numpy.where(x > 0.5)
u0[mask] = 2.0

# Plot of initial condition
plt.figure(figsize=(6.0, 4.0))
plt.title('Initial conditions')
plt.xlabel('x')
plt.ylabel('u')
plt.grid()
plt.plot(x, u0, color='C0', linestyle='-', linewidth=2)
plt.xlim(0.0, L)
plt.ylim(0.0, 2.1)

# Setting the coefficent matrix operator A
def lhs_operator(N, r):
    # Setup the diagonal of the operator.
    D = numpy.diag(1.0 * numpy.ones(N))
    # Setup the Neumann condition for the last element.
    D[-1, -1] = 1.0 + r / 2.0  # r = Courant number
    # upper diagonal of the operator.
    U = numpy.diag((r / 2) * numpy.ones(N - 1), k=1)
    # lower diagonal of the operator.
    L = numpy.diag((-1.0 * r / 2.0) * numpy.ones(N - 1), k=-1)
    # the final coeff matrix A
    A = D + U + L
    return A

# Setting the LHS vector with boundary conditions
def rhs_vector(u, r, qdx):
    b = u[1:-1]
    # Set Dirichlet condition.
    b[0] += u[0] * ( r / 2.0)
    # Set Neumann condition.
    b[-1] += qdx
    return b

# Solving the linear system of equation
def btcs_implicit(u0, nt, dt, dx, C, q):
    r = C * dt / dx
    # Create the implicit operator of the system.
    A = lhs_operator(len(u0) - 2, r)
    # Integrate in time.
    u = u0.copy()
    for n in range(nt):
        # Generate the right-hand side of the system.
        b = rhs_vector(u, r, q * dx)
        # Solve the system with scipy.linalg.solve.
        u[1:-1] = linalg.solve(A, b)
        # Apply the Neumann boundary condition.
        u[-1] = u[-2] + q * dx
    return u

# Specify the wave speed C, courant number r and q
r = 2
C = 1.0
q = 0
dt = r * dx / C # time-step size based of r, C and dx

# Compute velocity profile of the wave along the domain.
u = btcs_implicit(u0, nt, dt, dx, C, q)

# Plot of wave profile along the domain at final time step.
plt.figure(figsize=(6.0, 4.0))
plt.title('Wave profile at time step-10')
plt.xlabel('x')
plt.ylabel('u')
plt.grid()
plt.plot(x, u, color='C1', linestyle='-', linewidth=2)
plt.xlim(0.0, L)
plt.ylim(0.0, 2);

# calculating the analytical solution, u_an
t10 = 1 # time in seconds at tenth time step
xa = x-(C * t10)  # u_an = u0(x-Ct)
u_an = []
for i in range(len(x)):
  if x[i]-1<=0.5:
    u_an.append(1)
  else:
    u_an.append(2)

# Plot of analytical solution
plt.figure(figsize=(6.0, 4.0))
plt.title('analytical solution')
plt.xlabel('x')
plt.ylabel('u_an')
plt.grid()
plt.plot(x, u_an, color='C4', linestyle='-', linewidth=2)
plt.xlim(0.0, L)
plt.ylim(0, 2.1)

# Plot of L1 norm
norm = u_an - u
plt.figure(figsize=(6.0, 4.0))
plt.title('L1 norm')
plt.xlabel('x')
plt.ylabel('L1')
plt.grid()
plt.plot(x, norm, color='C2', linestyle='-', linewidth=2)
plt.xlim(0.0, L)
plt.ylim(0, 0.5)
