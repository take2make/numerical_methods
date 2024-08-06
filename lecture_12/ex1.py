from numpy import zeros, linspace
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import style, figure, axes

# y' = -y*sin(t) + exp(-cos(t));
# y(0) = 0;

def f(u, t):
    u_clipped = np.clip(u, -1e10, 1e10)
    f = -u_clipped*np.sin(t) + np.exp(-np.cos(t))
    return f

def euler_method(t_0, T, u_0, M):
    u = zeros(M + 1)
    u[0] = u_0
    tau = (T - t_0) / M
    t = linspace(t_0,T,M + 1)
    for m in range(M):
        u[m + 1] = u[m] + tau * f(u[m], t[m])
    return t, u

def erk2_method(t_0, T, u_0, M):
    u = zeros(M + 1)
    u[0] = u_0
    tau = (T - t_0) / M
    t = linspace(t_0,T,M + 1)
    for m in range(M):
        w_1 = f(u[m], t[m])
        w_2 = f(u[m] + tau * 2/3 * w_1, t[m] + tau * 2/3)
        u[m + 1] = u[m] + tau * (1/4 * w_1 + 3/4 * w_2)
    return t, u
        
t_0 = -1.; u_0 = 1; T = 50
M = 300

def main():
    #t, u = euler_method(t_0, T, u_0, M)
    t, u = erk2_method(t_0, T, u_0, M)
    plot(t, u)

def plot(t, u):
    style.use('dark_background')
    fig = figure()
    ax = axes(xlim=(t[0], t[-1]), ylim=(min(u)-3, max(u)+3))
    ax.plot(t, u, '-w', markersize=5)
    ax.set_title('Euler Method')
    plt.show()
    
if __name__ == "__main__":
    main()