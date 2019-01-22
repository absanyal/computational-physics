# Tue Jan 15 10:56:36 IST 2019

import numpy as np
import matplotlib.pyplot as plt

np.random.seed(1706)

sites = 100
J = 1
h = 0
T = 10
steps = 1000

upspin = 1
downspin = -1

#lattice = np.random.choice([upspin, downspin], size = sites)
lattice = np.random.choice([upspin, downspin], size=sites)
#lattice = np.zeros(sites, dtype = int)

print(lattice)


def energy_change(i):
    s = lattice[i]
    sm1 = lattice[i-1]
    sp1 = lattice[(i+1) % sites]
    e_i = -J*(sm1*s + s*sp1) - h*s

    s = -lattice[i]
    sm1 = lattice[i-1]
    sp1 = lattice[(i+1) % sites]
    e_f = -J*(sm1*s + s*sp1) - h*s

    return (e_f - e_i)


energy = 0
for i in range(sites):
    s = lattice[i]
    sp1 = lattice[(i+1) % sites]
    energy += -J*(s*sp1) - h * s
# print(energy)

energylist = [energy]
timelist = [0]
magnetization = [sum(lattice)/sites]

for t in range(1, steps):
    for i in range(sites):
        delta_e = energy_change(i)

        if (delta_e <= 0):
            lattice[i] *= -1
            energy += delta_e
        else:
            p = np.exp(-(1/T) * delta_e)
            r = np.random.uniform(0, 1)
            if (r <= p):
                lattice[i] *= -1
                energy += delta_e

    timelist.append(t)
    energylist.append(energy)
    magnetization.append(sum(lattice)/sites)
#    print(lattice)

    if (t > 0 and t % (steps/10) == 0):
        print("Done for", t, "steps.")

# print(energy)

print(lattice)
print("Initial magnetization", magnetization[0])
print("Average magnetization after thermalization",
      sum(magnetization[int(0.5 * steps):]) / (steps - int(0.5 * steps)))
plt.plot(timelist, energylist)
plt.show()
plt.ylim(-1, 1)
plt.plot(timelist, magnetization)
