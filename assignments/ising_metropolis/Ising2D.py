# Tue Jan 15 10:56:42 IST 2019

import numpy as np
import matplotlib.pyplot as plt

np.random.seed(1706)

J = 1
h = 0
T = 2.3
steps = 1000

length = 20
breadth = 20

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

energylist = [energy]
timelist = [0]
magnetization = [sum(sum(lattice))/(length*breadth)]

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

    timelist.append(t)
    energylist.append(energy)
    magnetization.append(sum(sum(lattice))/(length*breadth))
#    print(lattice)

    if (t > 0 and t % (steps/10) == 0):
        print("Done for", t, "steps.")

# print(energy)

print(lattice)
print("Initial magnetization", magnetization[0])
print("Average magnetization after thermalization",
      sum(magnetization[int(0.5 * steps):]) / (steps - int(0.5 * steps)))

plt.ylim(-1, 1)
plt.plot(timelist, energylist)
plt.show()
plt.plot(timelist, magnetization)
plt.show()
