from numpy import zeros, linspace
import matplotlib.pyplot as plt
from matplotlib.pyplot import style, figure, axes

def f(u, t, m_sun, G):
    """
    Function f prepares an array containing the elements of the vector function
    that defines the right-hand side of the system of ODEs
    """
    f = zeros(4)
    f[0] = u[2]
    f[1] = u[3]
    f[2] = -G * m_sun * u[0] / (u[0]**2 + u[1]**2)**(3/2)
    f[3] = -G * m_sun * u[1] / (u[0]**2 + u[1]**2)**(3/2)
    return f

def solve_ode(t_0, T, x_0, y_0, v_x_0, v_y_0, G, m_sun, M, s, a, b, c):
    """
    Solve the system of ODEs using the ERK4 scheme
    """
    tau = (T - t_0) / M
    t = linspace(t_0, T, M + 1)
    u = zeros((M + 1, 4))
    u[0] = [x_0, y_0, v_x_0, v_y_0]

    for m in range(M):
        w = zeros((s, 4))
        for k in range(s):
            adjustment_1 = zeros(4)
            for l in range(k):
                adjustment_1 = adjustment_1 + a[k, l] * w[l]
            w[k] = f(u[m] + tau * adjustment_1, t[m] + tau * c[k], m_sun, G)
        adjustment_2 = zeros(4)
        for k in range(s):
            adjustment_2 = adjustment_2 + b[k] * w[k]
        u[m + 1] = u[m] + tau * adjustment_2

    return t, u

def plot_trajectory(t, u):
    """
    Plot the trajectory of Earth
    """
    style.use('dark_background')
    fig = figure()
    ax = axes(xlim=(-2 * 10**11, 2 * 10**11), ylim=(-2 * 10**11, 2 * 10**11))
    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.plot(0, 0, 'yo', markersize=15)
    ax.plot(u[:, 0], u[:, 1], '-w', markersize=5)
    ax.plot(u[-1, 0], u[-1, 1], color='w', marker='o', markersize=7)
    ax.set_title("Earth's Trajectory")
    plt.show()

# Define the input data for the problem
t_0 = 0.
T = 365.25 * 24 * 60 * 60
x_0 = 147098291 * 10**3
y_0 = 0.
v_x_0 = 0.
v_y_0 = 30.4 * 10**3
G = 6.674301515151515 * 10**(-11)
m_sun = 1.98847 * 10**30

# Define the number of intervals in the grid
M = 365

# Define the ERK scheme to be used for calculations
# s = 2  # ERK2
# b = zeros(s)
# a = zeros((s, s))
# c = zeros(s)
# b[0] = 1/4
# b[1] = 3/4
# a[1, 0] = 2/3
# c[0] = 0.
# c[1] = 2/3

s = 4 # ERK4
b = zeros(s); a = zeros((s,s)); c = zeros(s)
b[0] = 1/6; b[1] = 1/3; b[2] = 1/3; b[3] = 1/6
a[1,0] = 1/2; a[2,0] = 0.; a[2,1] = 1/2; a[3,0] = 0.; a[3,1] = 0.; a[3,2] = 1.
c[0] = 0.; c[1] = 1/2; c[2] = 1/2; c[3] = 1.

def main():
    # Solve the system of ODEs and plot the trajectory
    t, u = solve_ode(t_0, T, x_0, y_0, v_x_0, v_y_0, G, m_sun, M, s, a, b, c)
    plot_trajectory(t, u)
    
if __name__ == '__main__':
    main()