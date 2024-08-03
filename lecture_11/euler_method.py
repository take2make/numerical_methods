from numpy import zeros, linspace
import matplotlib.pyplot as plt
from matplotlib.pyplot import style, figure, axes
from celluloid import Camera

# This code solves Newton's equations for Earth's movement 
# around the Sun using the Euler method.

# The function f prepares an array containing the elements of the vector function,
# which defines the right-hand side of the system of ordinary differential equations being solved
def f(u,t,m_sun,G):
    f = zeros(4)
    f[0] = u[2]
    f[1] = u[3]
    f[2] = - G*m_sun*u[0]/(u[0]**2 + u[1]**2)**(3/2)
    f[3] = - G*m_sun*u[1]/(u[0]**2 + u[1]**2)**(3/2)
    return f

# Definition of the problem's input data
t_0 = 0.; T = 365.25*24*60*60
x_0 = 147098291*10**3; y_0 = 0.
v_x_0 = 0.; v_y_0 = 30.4*10**3
G = 6.674301515151515*10**(-11)
m_sun = 1.98847*10**30

# Definition of the number of intervals in the grid
M = 3650

# Definition of the grid
tau = (T - t_0)/M
t = linspace(t_0,T,M + 1)
# Allocating memory for the array of grid values of the solution to the system of ODEs
# The grid values of the solution corresponding to the time t_m are stored in the row with number m of this array
u = zeros((M + 1,4))

# Setting initial conditions
# (written in row number 0 of the u array)
u[0] = [x_0, y_0, v_x_0, v_y_0]

# Implementation of the Euler method
for m in range(M):
    u[m + 1] = u[m] + tau*f(u[m],t[m],m_sun,G)

# Plotting the solution
def paint():
    style.use('dark_background')
    fig = figure()
    ax = axes(xlim=(-2*10**11,2*10**11), ylim=(-2*10**11,2*10**11))
    ax.set_aspect('equal'); ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.plot(0,0,'yo',markersize=15)
    ax.plot(u[:,0],u[:,1],'-w',markersize=5)
    ax.plot(u[M,0],u[M,1], color='w', marker='o', markersize=7)
    ax.set_title('Earth\'s Trajectory')
    plt.show()
    
# Animation of the solution
def animate():
    style.use('dark_background')
    fig = figure()
    camera = Camera(fig)
    ax = axes(xlim=(-2*10**11,2*10**11), ylim=(-2*10**11,2*10**11))
    ax.set_aspect('equal'); ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_title('Earth\'s Trajectory')
    for m in range(M + 1):
        # Plotting the static Sun at the origin
        ax.plot(0,0,'yo',markersize=15)
        # Plotting the Earth at time t[m]
        ax.plot(u[m,0],u[m,1], color='w', marker='o', markersize=7)
        # Plotting the path traveled by the Earth up to time t[m]
        ax.plot(u[:m+1,0], u[:m+1,1], color='w', ls='--', lw=2)
        camera.snap()
    animation = camera.animate(interval=15, repeat=False, blit=True)
    plt.show()

def main():
    paint()

if __name__ == "__main__":
    main()