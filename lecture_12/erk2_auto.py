from numpy import zeros, linspace
import matplotlib.pyplot as plt
from matplotlib.pyplot import style, figure, axes

# The function f prepares an array containing the elements of the vector function
# that defines the right-hand side of the system of ODEs being solved
def f(u, lambd):
    f = zeros(2)
    f[0] = lambd * u[0] * (u[1] - u[0])
    f[1] = 1
    return f

# Definition of the input data for the problem
t_0 = -1.; T = 2.
u_0 = 3.; lambd = 10.

# Definition of the number of intervals in the grid
# on which the approximate solution will be sought
M = 250

# Definition of the grid
tau = (T - t_0) / M
t = linspace(t_0, T, M + 1)

# Memory allocation for the array of grid values of the solution to the system of ODEs
# The grid values corresponding to the time moment t_m are stored in the row with number m of this array
u = zeros((M + 1, 2))

# Setting the initial conditions
# (written to the row with number 0 of the u array)
u[0] = [u_0, t_0]

# # Implementation of the Euler scheme
# for m in range(M):
#     u[m + 1] = u[m] + tau * f(u[m], lambd)

# Implementation of the ERK2 scheme
for m in range(M):
    w_1 = f(u[m], lambd)
    w_2 = f(u[m] + tau * 2 / 3 * w_1, lambd)
    u[m + 1] = u[m] + tau * (1 / 4 * w_1 + 3 / 4 * w_2)

# Plotting the solution
def paint():
    style.use('dark_background')
    ax = axes(xlim=(-1, 2), ylim=(0, 3))
    ax.set_aspect('equal')
    ax.set_xlabel('t')
    ax.set_ylabel('u')
    ax.plot(u[:, 1], u[:, 0], '-y', lw=3)
    ax.set_title('Graph of u(t)')
    plt.show()

def main():
    paint()
    
if __name__ == '__main__':
    main()

# Listing of the program that implements the solution to the autonomous system of ODEs
# after the autonomousization procedure has been performed
