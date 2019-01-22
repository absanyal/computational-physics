# Tue Jan 15 10:56:42 IST 2019

import numpy as np
import matplotlib.pyplot as plt

#np.random.seed(617)

J = 1
h = 0
# T = 2.3
steps = 1000

length = 10
breadth = 10

upspin = 1
downspin = -1

#lattice = np.zeros((length, breadth), dtype = int)
lattice = np.random.choice([upspin, downspin], size=(length, breadth))
print(lattice)


def energy_change(i, j):
    s = lattice[i, j]
    sym1 = lattice[i-1, j]
    syp1 = lattice[(i+1) % length, j]
    sxm1 = lattice[i, j - 1]
    sxp1 = lattice[i, (j+1) % breadth]
    e_i = -J*(s * (sxp1 + sxm1 + syp1 + sym1)) - h * s

    s = - lattice[i, j]
    sym1 = lattice[i-1, j]
    syp1 = lattice[(i+1) % length, j]
    sxm1 = lattice[i, j - 1]
    sxp1 = lattice[i, (j+1) % breadth]
    e_f = -J*(s * (sxp1 + sxm1 + syp1 + sym1)) - h * s

    return (e_f - e_i)

energy = 0
for i in range(length):
    for j in range(breadth):
        s = lattice[i, j]
        syp1 = lattice[(i+1) % length, j]
        sxp1 = lattice[i, (j+1) % breadth]
        energy += -J*(s * (sxp1 + syp1)) - h * s

# print(energy)


timelist = [0]

avg_m_list = []
avg_e_list = []

T = 4
T_list = []
while (T >= 0.1):
    magnetization = [sum(sum(lattice))/(length*breadth)]
    energy = 0
    for i in range(length):
        for j in range(breadth):
            s = lattice[i, j]
            syp1 = lattice[(i+1) % length, j]
            sxp1 = lattice[i, (j+1) % breadth]
            energy += -J*(s * (sxp1 + syp1)) - h * s
    energylist = [energy]
    for t in range(1, steps):
        for i in range(length):
            for j in range(breadth):
                delta_e = energy_change(i, j)

                if (delta_e <= 0):
                    lattice[i, j] *= -1
                    energy += delta_e
                else:
                    p = np.exp(-(1/T) * delta_e)
                    r = np.random.uniform(0, 1)
                    if (r <= p):
                        lattice[i, j] *= -1
                        energy += delta_e

        # timelist.append(t)
        energylist.append(energy)
        magnetization.append(sum(sum(lattice))/(length*breadth))
    #    print(lattice)

        # if (t > 0 and t % (steps/10) == 0):
        #     print("Done for", t, "steps.")
    avg_m = sum(magnetization[int(0.5 * steps):])/(steps - int(0.5 * steps))
    avg_e = sum(energylist[int(0.5 * steps):])/(steps - int(0.5 * steps)) \
            / (length * breadth)
    avg_m_list.append(avg_m)
    avg_e_list.append(avg_e)
    T_list.append(T)
    print(round(T, 2), '\t', round(avg_m, 5), '\t', round(avg_e, 5))
    T -= 0.1


plt.plot(T_list, avg_e_list)
plt.title('Energy vs T')
plt.ylabel('Average energy per site')
plt.xlabel('Temperature')
plt.savefig('E_vs_T.pdf')
plt.show()
plt.ylim(-1, 1)
plt.plot(T_list, avg_m_list)
plt.title('Magnetization vs T')
plt.ylabel('Average magnetization per site')
plt.xlabel('Temperature')
plt.savefig('m_vs_T.pdf')
