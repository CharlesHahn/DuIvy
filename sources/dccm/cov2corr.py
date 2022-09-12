import pandas as pd
import numpy as np
import math 
import matplotlib
from matplotlib.pylab import *
import matplotlib.pyplot as plt

lis = ["#80FFFF","#A8FFFF","#D4FFFF","#FFFCFF", "#FFFCFF","#FFD4FF", "#FFA8FF", "#FF80FF"]
# lis = ["#FF80FF", "#FFA8FF", "#FFD4FF", "#FFFCFF", "#FFFCFF", "#D4FFFF", "#A8FFFF", "#80FFFF"]
lis.reverse()
cm = matplotlib.colors.ListedColormap(lis)

#Load the output file of GROMACS covar ascii.
covar = pd.read_csv('hhh.dat', sep=' ', header=None)

#Parsing the covariance file
resnum = int(math.sqrt((len(covar.index))/3))

all_results = pd.DataFrame()

for i in range(0,resnum):
    three_step = pd.DataFrame()
    print((i*resnum)*3,int(len(covar)/resnum)*(i+1),resnum)
    for j in range((i*resnum)*3,int(len(covar)/resnum)*(i+1),resnum):
        df = covar[j:resnum+j]
        df2 = df.reset_index(drop=True)
        three_step = pd.concat([three_step,df2], ignore_index=True, axis=1)
    
    all_results = pd.concat([all_results,three_step], ignore_index=True, axis=0)

print(all_results.shape)
all_results['sum'] = all_results.sum(axis=1)
print(all_results.shape)
a=all_results['sum'].to_numpy()
print(a)
cov_matrix = a.reshape(resnum,resnum)
print(cov_matrix)

#Convert the covariance matrix to cross-correlation
corr=np.zeros((resnum,resnum))
for i in range(0,resnum):
    for j in range(0,resnum):
        corr[i,j] = cov_matrix[i,j]/math.sqrt(cov_matrix[i,i]*cov_matrix[j,j])


#Save the cross-correlation matrix as csv file for further usage
np.savetxt("corr.csv", corr, delimiter=" ", fmt='%s')

#Draw the graph
# file=np.loadtxt('corr.csv')
# data= np.array(file)
fig, ax = plt.subplots(figsize=(11,9))
data = pd.read_csv("coor2corr.csv")
print(data)


#Set the vmin and vmax values according to your interests
#You can change the cmap styles (PiYG, PRGn, BrBG, PuOr, RdGy, RdBu,
#RdYlBu,RdYlGn, Spectral, coolwarm, bwr, seismic, twilight.. )
"""
im = ax.imshow(data, cmap=cm, vmin=-1, vmax=1, origin='lower')
cbar = ax.figure.colorbar(im, ax=ax)
cbar.ax.tick_params(labelsize=24)
for tick in cbar.ax.get_yticklabels():
    tick.set_fontname("Times New Roman")
ax.set_title("Cross Correlation", fontname= "Times New Roman", fontsize=42, pad=15)
ax.tick_params(labelsize=24)
for tick in ax.get_xticklabels():
    tick.set_fontname("Times New Roman")
for tick in ax.get_yticklabels():
    tick.set_fontname("Times New Roman")
fig.tight_layout()
plt.show()

plt.clf()
im = plt.imshow(data,cmap=cm, vmin=-1, vmax=1, origin="lower", interpolation="gaussian")
cb = plt.colorbar(im)
plt.tight_layout()
plt.show()
"""

plt.clf()
im = plt.contourf(data, cmap=cm)
plt.contour(data, linewidths=0.2, levels=[-1.00, -0.75, -0.50, -0.25, 0.25, 0.50, 0.75, 1.00])
cb = plt.colorbar(im)
plt.xlabel("Residue No.")
plt.ylabel("Residue No.")
plt.title("Residue Cross Correlation")
plt.tight_layout()
plt.show()





#Save the graph
# fig.savefig('corr.png', dpi=300)

 
