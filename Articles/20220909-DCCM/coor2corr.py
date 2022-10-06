## author : charlie
## date : 20220831

import numpy as np
import pandas as pd
import math
from matplotlib import pyplot as plt
from matplotlib import colors

lis = ["#80FFFF","#A8FFFF","#D4FFFF","#FFFCFF", "#FFFCFF","#FFD4FF", "#FFA8FF", "#FF80FF"]
lis.reverse()
cm = colors.ListedColormap(lis)

with open("alpha_c.gro", 'r') as fo:
    lines = fo.readlines()
atom_number = int(lines[1].strip())
coordinates = []

for i in range(0, len(lines), atom_number+3):
    coordinates.append([])
    for j in range(i+2,i+2+atom_number):
        x = float(lines[j][20:28])*10
        y = float(lines[j][28:36])*10
        z = float(lines[j][36:44])*10
        coordinates[-1].append([x, y, z])
xyz = np.array(coordinates)
print(xyz.shape)


average = []
with open("alpha_c_ref.gro", 'r') as fo:
    lines = fo.readlines()[2:-1]
for i in range(atom_number):
    x = float(lines[i][20:28])*10
    y = float(lines[i][28:36])*10
    z = float(lines[i][36:44])*10
    average.append([x, y, z])

time_number = len(coordinates)

# kabschalgorithm 点云算法修正坐标 fit
# fitxyz
allxyz = np.array(coordinates)
for n in range(allxyz.shape[0]):
    cercoord_A=np.mean(np.matrix(allxyz[0]),axis=0)
    cercoord_B=np.mean(np.matrix(allxyz[n]),axis=0)
    A_c=allxyz[0]-cercoord_A
    B_c=allxyz[n]-cercoord_B
    H=B_c.T*A_c
    U,S,VT=np.linalg.svd(H)
    d=np.sign(np.linalg.det(VT.T*U.T))
    diag=np.diag([1,1,d])
    R=VT.T * diag *  U.T
    # 坐标以A为主
    T=-R*cercoord_B.T+cercoord_A.T
    # 修正之后的B坐标
    B_calc=(R*allxyz[n].T+T).T
    allxyz[n]=B_calc

# 位置偏移量是基于平均结构计算的，而不是参考结构
atomcoord_mean=[]
for i in range(allxyz.shape[1]):
    x=[]
    y=[]
    z=[]
    for j in range(allxyz.shape[0]):
        x.append(allxyz[j][i][0])
        y.append(allxyz[j][i][1])
        z.append(allxyz[j][i][2])
    atomcoord_mean.append([np.mean(x),np.mean(y),np.mean(z)])
delta0=np.array(atomcoord_mean)   
offset = allxyz - delta0
print("fit done")


covariance = np.zeros((atom_number, atom_number))
for i in range(atom_number):
    for j in range(atom_number):
        covariance[i, j] = np.mean([np.dot(offset[t][i], offset[t][j]) for t in range(time_number)])

# covariance = np.array(covariance)
# cmax, cmin = covariance.max(), covariance.min()
# covariance = (covariance - cmin)/(cmax-cmin)

print(covariance)
plt.imshow(covariance, cmap=cm, origin="lower")
plt.show()

corr=np.zeros((atom_number,atom_number))
for i in range(0,atom_number):
    for j in range(0,atom_number):
        corr[i,j] = covariance[i,j]/np.sqrt(covariance[i,i]*covariance[j,j])
print(corr)

fig, ax = plt.subplots(figsize=(11,9))
im = ax.imshow(corr, cmap=cm, vmin=-1, vmax=1, origin='lower')
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

df = pd.DataFrame(corr)
df.to_csv("coor2corr.txt", sep=" ")



