from numpy import zeros, linspace, log, sqrt, sum
from matplotlib.pyplot import style, figure, axes

import matplotlib.pyplot as plt

# The function f prepares an array containing the elements of the vector function
# that defines the right-hand side of the system of ODEs being solved
def f(u, m_sun, G):
    f = zeros(4)
    f[0] = u[2]
    f[1] = u[3]
    f[2] = - G * m_sun * u[0] / (u[0]**2 + u[1]**2)**(3/2)
    f[3] = - G * m_sun * u[1] / (u[0]**2 + u[1]**2)**(3/2)
    return f

# The function implements the solution of the system of ODEs
# on a grid with M_0*r**s intervals using the ERK2 scheme
def ODESolving(t_0, T, x_0, y_0, v_x_0, v_y_0, M_0, s, r, m_sun, G):
    # Input parameters:
    # t_0, T - initial and final time moments
    # x_0, y_0, v_x_0, v_y_0 - initial conditions
    # M_0 - number of intervals in the base grid in time
    # s - grid number on which the solution is computed
    # (if s = 0, then the solution is computed on the base grid)
    # r - grid refinement coefficient
    # m_sun, G - problem parameters

    # Output parameter:
    # u_basic - array containing grid values
    # of the solution of the ODE system only at nodes
    # coinciding with nodes of the base grid

    # Formation of a refined grid with a refinement factor of r^s:

    # Calculation of the number of intervals on the grid with number s
    M = M_0 * r**s
    # Determination of the step of the refined grid
    tau = (T - t_0) / M
    # Determination of the refined grid
    t = linspace(t_0, T, M + 1)

    # Memory allocation for the array of grid values of the ODE solution,
    # which will store grid values from nodes
    # coinciding with nodes of the base grid
    u_basic = zeros((M_0 + 1, 4))

    # Memory allocation for the array of grid values
    # of the solution on the refined grid
    # The m-th row of this array stores grid values of the solution
    # corresponding to the time moment t_m
    u = zeros((M + 1, 4))

    # Setting the initial condition
    u[0] = [x_0, y_0, v_x_0, v_y_0]

    # Implementation of the ERK2 scheme
    for m in range(M):
        w_1 = f(u[m], m_sun, G)
        w_2 = f(u[m] + tau * 2/3 * w_1, m_sun, G)
        u[m + 1] = u[m] + tau * (1/4 * w_1 + 3/4 * w_2)

    # Selecting grid values from nodes
    # coinciding with nodes of the base grid from the u array
    for m in range(M_0 + 1):
        u_basic[m] = u[m * r**s]

    return u_basic

# Definition of the input data for the problem
t_0 = 0.; T = 365.25 * 24 * 60 * 60
x_0 = 147098291 * 10**3; y_0 = 0.
v_x_0 = 0.; v_y_0 = 30.4 * 10**3
G = 6.674301515151515 * 10**(-11)
m_sun = 1.98847 * 10**30

# Definition of the number of intervals in the BASE grid
# on which the approximate solution will be sought
M = 36

# Number of grids on which the approximate solution is sought
S = 7
# Grid refinement coefficient
r = 2
# Theoretical parameters of the scheme
p = 2; q = 1

# Memory allocation for arrays of grid values
# of the ODE solutions on different grids with numbers s = 0,...,S-1,
# which store grid values of the solution from nodes
# coinciding with nodes of the base grid
U = zeros((S, S, M + 1, 4))

# "Big loop" that recalculates the solution S times
# on a sequence of refining grids
# The array of grid values of the solution contains only
# grid values from nodes coinciding with nodes of the base grid
for s in range(S):
    U[s, 0, :, :] = ODESolving(t_0, T, x_0, y_0, v_x_0, v_y_0, M, s, r, m_sun, G)

# Memory allocation for arrays of errors R,
# relative errors R_rel, and effective orders of accuracy p_eff
R = zeros((S, S, M + 1, 4))
R_rel = zeros((S, S))
p_eff = zeros((S, S))

for s in range(1, S):
    for l in range(s):
        R[s, l, :, :] = (U[s, l, :, :] - U[s-1, l, :, :]) / (r**(p + l*q) - 1)
        U[s, l+1, :, :] = U[s, l, :, :] + R[s, l, :, :]
        R_rel[s, l] = sqrt(sum(R[s, l, :, :]**2)) / sqrt(sum(U[s, l+1, :, :]**2)) * 100

for s in range(2, S):
    for l in range(s-1):
        p_eff[s, l] = log(sqrt(sum(R[s-1, l, :, :]**2)) / sqrt(sum(R[s, l, :, :]**2))) / log(r)

# The function prints a formatted table
def PrintTriangular(A, i):
    print(' ', end=' ')
    for l in range(len(A)):
        print(' p={0:<4d}'.format(p + l*q), end=' ')
    print()
    for m in range(len(A)):
        print('s={0:<2d}'.format(m), end=' ')
        for l in range(m + 1 - i):
            print('{0:5.2f}'.format(A[m, l]), end=' ')
        print()
    print()

print('Table of relative error estimates (in percentages):')
PrintTriangular(R_rel, 1)
print('Table of effective orders of accuracy:')
PrintTriangular(p_eff, 2)

# Plotting the solution obtained on the grid with number S-1
# (only nodes coinciding with nodes of the base grid are marked)
style.use('dark_background')

fig1 = figure()
ax1 = axes(xlim=(-2*10**11, 2*10**11), ylim=(-2*10**11, 2*10**11))
ax1.set_aspect('equal'); ax1.set_xlabel('x'); ax1.set_ylabel('y');
ax1.plot(0, 0, 'yo', markersize=15)
ax1.plot(U[S-1, 1, :, 0], U[S-1, 1, :, 1], '-w', markersize=5)
ax1.plot(U[S-1, 1, M, 0], U[S-1, 1, M, 1], color='w', marker='o', markersize=7)
ax1.set_title('Earth Trajectory')
plt.show()

# Plotting the dependence of the error on the number of grid intervals
fig2 = figure()
ax2 = axes()
ax2.plot([r**s*M for s in range(1, S)], [sqrt(sum(R[s, 0, :, :]**2)) for s in range(1, S)], '-wo')
ax2.set_xscale('log'); ax2.set_yscale('log')
plt.show()
