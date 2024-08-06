from numpy import zeros, linspace
import matplotlib.pyplot as plt
from matplotlib.pyplot import style, figure, axes

# The function f prepares an array containing the elements of the vector function
# that defines the right-hand side of the system of ODEs being solved
def f(u, t, m_sun, G):
    f = zeros(4)
    f[0] = u[2]
    f[1] = u[3]
    f[2] = - G * m_sun * u[0] / (u[0]**2 + u[1]**2)**(3/2)
    f[3] = - G * m_sun * u[1] / (u[0]**2 + u[1]**2)**(3/2)
    return f

# Definition of the input data for the problem
t_0 = 0.; T = 365.25 * 24 * 60 * 60
x_0 = 147098291 * 10**3; y_0 = 0.
v_x_0 = 0.; v_y_0 = 30.4 * 10**3
G = 6.674301515151515 * 10**(-11)
m_sun = 1.98847 * 10**30

# Definition of the number of intervals in the grid
# on which the approximate solution will be sought
M = 365

# Definition of the grid
tau = (T - t_0) / M
t = linspace(t_0, T, M + 1)

# Memory allocation for the array of grid values of the solution to the system of ODEs
# The m-th row of this array stores the grid values of the solution
# corresponding to the time moment t_m
u = zeros((M + 1, 4))

# Setting the initial conditions
# (written in the 0-th row of the u array)
u[0] = [x_0, y_0, v_x_0, v_y_0]

# Implementation of the ERK2 scheme
def erk2():
    for m in range(M):
        w_1 = f(u[m], t[m], m_sun, G)
        w_2 = f(u[m] + tau * 2/3 * w_1, t[m] + tau * 2/3, m_sun, G)
        u[m + 1] = u[m] + tau * (1/4 * w_1 + 3/4 * w_2)

# Plotting the solution 
def paint():
    style.use('dark_background')
    fig = figure()
    ax = axes(xlim=(-2 * 10**11, 2 * 10**11), ylim=(-2 * 10**11, 2 * 10**11))
    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.plot(0, 0, 'yo', markersize=15)
    ax.plot(u[:, 0], u[:, 1], '-w', markersize=5)
    ax.plot(u[M, 0], u[M, 1], color='w', marker='o', markersize=7)
    ax.set_title('Earth\'s Trajectory')
    plt.show()
    
def main():
    erk2()
    paint()
    
if __name__ == '__main__':
    main()