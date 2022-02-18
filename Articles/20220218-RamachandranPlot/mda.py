import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran

u = mda.Universe("pro.gro", "test.xtc")
r = u.select_atoms("protein")
R = Ramachandran(r).run()
fig, ax = plt.subplots()
R.plot(ax=ax, color='k', marker='s', ref=True)
plt.show()


