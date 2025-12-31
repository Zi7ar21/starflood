from galpy import df
from galpy import potential
from galpy.actionAngle import actionAngleStaeckel
import matplotlib.pyplot as plt

print("This script is incomplete, feel free to make any contributions/suggestions on the Starflood repository: https://github.com/Zi7ar21/starflood")

print([p for p in dir(potential) if 'Potential' in p])

aAS = actionAngleStaeckel(pot=potential.MWPotential2014, delta=0.45, c=True)

qdfS = df.quasiisothermaldf(1.0/3.0, 0.2, 0.1, 1.0, 1.0, pot=potential.MWPotential2014, aA=aAS, cutcounter=True)

vs = qdfS.sampleV(1.0, 0.0, n=1000)

plt.hist(vs[:,1], density=True, histtype='step', bins=101, range=[0.0, 1.5])

plt.show()
