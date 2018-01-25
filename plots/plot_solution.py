import numpy as np
import matplotlib.pyplot as plt

# MKS
sec = 1.0
cm = 1e-2
km = 1e3
pc = 3.086e16
kpc = 1e3 * pc

def plot_solution(Q0, v, L, D):
    z = np.linspace(-L, L, 100)
    w = -v * L / D
    N = Q0 / 2. / v
    N *= 1. - np.exp(w * (1. - np.abs(z) / L))
    N /= 1. - np.exp(w)
    plt.plot(z / kpc, N, 'k:')

for i in range(0, 501, 100):
    print i
    filename = 'output/N_size_129_t_' + str(i) + '.txt'
    z, N = np.loadtxt(filename,skiprows=1,usecols=(0, 1),unpack=True)
    plt.plot(z, N, 'r')
    filename = 'output/N_size_257_t_' + str(i) + '.txt'
    z, N = np.loadtxt(filename,skiprows=1,usecols=(0, 1),unpack=True)
    plt.plot(z, N, 'b')


plot_solution(1e-30, 30 * km / sec, 4 * kpc, 3e28 * cm**2 / sec)

#plt.yscale('log')
plt.ylim([0, 2e-35])
plt.show()
