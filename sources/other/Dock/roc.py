# author: charlie 
# roc test 

import os
import matplotlib.pyplot as plt


def y():
	for i in range(1000):
		yield str(i)

p = 0
n = 0
tpr_li = []
fpr_li = []


roc_dic = {}
yie = y()
for filename in os.listdir():
	if '.log' in filename and filename[0] == 't':
		with open(filename, 'r') as fo:
			content = fo.read()
		lines = content.strip().split('\n')
		for i in range(-10, -1):
			energy = lines[i].strip().split()[1]
			roc_dic[next(yie) + 't'] = float(energy)
			p += 1

	elif '.log' in filename and filename[0] == 'f':
		with open(filename, 'r') as fo:
			content = fo.read()
		lines = content.strip().split('\n')
		for i in range(-10, -1):
			energy = lines[i].strip().split()[1]
			roc_dic[next(yie) + 'f'] = float(energy)
			n += 1

roc = sorted(roc_dic.items(), key = lambda x: x[1])

# print(p, n)

for n in range(len(roc)-1):
	gate = (roc[n][1] + roc[n+1][1] ) * 0.5
	tp = 0
	fp = 0
	fn = 0
	tn = 0
	# print(roc[n][1], gate)
	for ro in roc:
		if 't' in ro[0] and ro[1] <= gate:
			tp += 1
		elif 'f' in ro[0] and ro[1] <= gate:
			fp += 1
		elif 't' in ro[0] and ro[1] > gate:
			fn += 1
		elif 'f' in ro[0] and ro[1] > gate:
			tn += 1
	#print(tp, fp, fn, tn)
	tpr = tp/(tp + fn)
	fpr = fp/(fp + tn)
	# print(tpr, fpr)
	tpr_li.append(tpr)
	fpr_li.append(fpr)



plt.plot(fpr_li, tpr_li)
plt.plot([0.1*i for i in range(11)], [0.1*i for i in range(11)])
plt.legend(['tpr&fpr', 'y=x'])
plt.title('ROC of vina')
plt.xlabel('fpr')
plt.ylabel('tpr')
plt.show()

