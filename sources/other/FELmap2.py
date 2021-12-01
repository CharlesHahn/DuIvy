# author : charlie

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import matplotlib.mlab as mlab
from matplotlib import pylab as pylab
import sys
from scipy import interpolate


myparams = {
	'axes.labelsize': '12',
	'xtick.labelsize': '12',
	'ytick.labelsize': '12',
	'lines.linewidth': 1,
	'legend.fontsize': '12',
	'font.family': 'Times New Roman'
}
pylab.rcParams.update(myparams)




def getdata(filename):
	with open(filename, 'r') as fo:
		content = fo.read()

	x, y, z = [], [], []
	lines = content.strip('\n').split('\n')
	for line in lines:
		x.append( float(line.split()[0] ) )
		y.append( float(line.split()[1] ) )
		z.append( float(line.split()[2] ) )

	return x, y, z 



def plot3D(x, y, z, ip):
	x = sorted( list( set(x)) ) 
	y = sorted( list( set(y)) )
	x = np.array(x)
	y = np.array(y)
	z = np.array(z)
	# z = z.reshape( (len(x), len(y) ) )
	# print(z.shape)

	# 插值
	xnew = np.linspace(min(x), max(x), ip*len(x))
	ynew = np.linspace(min(y), max(y), ip*len(y))
	f = interpolate.interp2d(x, y, z, kind="cubic")
	znew = f(xnew, ynew)
	x, y= xnew, ynew
	z = []
	for lis in znew:
		z += list(lis)
	# print(len(z))
	z = np.array(z)

	z = z.reshape( (len(x), len(y) ) )
	x, y = np.meshgrid(x, y)

	fig = plt.figure()
	ax = fig.gca(projection='3d')

	surf = ax.plot_surface(x, y, z, alpha = 0.9, cmap = cm.coolwarm, linewidth = 0, antialiased = False  )
	cset = ax.contourf(x, y, z, zdir='z', offset= -10, cmap=cm.coolwarm)

	ax.zaxis.set_major_locator(LinearLocator(10))

	ax.zaxis.set_major_formatter(FormatStrFormatter('%.00f'))
	ax.yaxis.set_major_formatter(FormatStrFormatter('%.02f'))
	ax.xaxis.set_major_formatter(FormatStrFormatter('%.02f'))
	ax.set_zlim(-10, )
	ax.set_xlabel('Gyrate (nm)')
	ax.set_ylabel('RMSD (nm)')
	ax.set_zlabel('Free Energy')
	fig.colorbar(surf, shrink=0.5, aspect=5)

	plt.show()



def plot2d(x, y, z):

	x = sorted( list( set(x)) ) 
	y = sorted( list( set(y)) )
	x = np.array(x)
	y = np.array(y)
	z = np.array(z)
	z = z.reshape( (len(x), len(y) ) )
	x, y = np.meshgrid(x, y)

	
	fig, ax = plt.subplots()

	norm = cm.colors.Normalize(vmax=20, vmin=-10 )
	cset = ax.contour(x, y, z, 40, norm = norm )

	ax.set_xlabel('gyrate')
	ax.set_ylabel('rmsd')
	ax.yaxis.set_major_formatter(FormatStrFormatter('%.02f'))
	ax.xaxis.set_major_formatter(FormatStrFormatter('%.02f'))

	plt.show()



def main():
	info = """
- python FELmap2.py file2draw 
		(-d2 : 2d plot) 
		(-ip+num : interpolate to num times more data)
	eg. FELmap2.py filename -ip10 
	"""
	print(info)
	plot_style = 3
	ip = 1 # 不插值
	for arg in sys.argv[1:]:
		if '-d' == arg[:2] and '2' == arg[2]:
			plot_style = 2
		elif '-ip' == arg[:3]:
			ip = int(arg[3:])
		else:
			filename = arg

	x, y, z = getdata(filename)
	if plot_style == 3:
		plot3D(x, y, z, ip)
	elif plot_style == 2:
		plot2d(x, y, z)



if __name__ == '__main__':
	main()
