import matplotlib.pyplot as plt
import numpy as np

# Exemple de courbe de rotation (données fictives pour illustration)
r = np.linspace(0.5, 10, 50)
vobs = 80 * np.sqrt(r) / np.sqrt(r + 1) + 20
vmodel = vobs * 0.95 + np.random.normal(0, 3, len(r))

plt.figure(figsize=(8, 5))
plt.plot(r, vobs, 'o', label='Vobs')
plt.plot(r, vmodel, '-', label='Vmodel')
plt.xlabel('R (kpc)')
plt.ylabel('V (km/s)')
plt.title('Exemple de courbe de rotation')
plt.legend()
plt.grid(True)
plt.savefig('results/example_rotation_curve.png')
print('Plot sauvegardé dans results/example_rotation_curve.png')
