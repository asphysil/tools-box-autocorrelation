import numpy as np 

from matplotlib import pyplot as plt 



#data = np.genfromtxt('md-vac.dat')

data = np.genfromtxt('md-ph-dos.dat')
#data = np.genfromtxt('expanded_data.txt')
#data = np.genfromtxt('test-ref.dat')
x=data[:,0]
print(len(x), len(x)/2 )
y=data[:,1]

plt.plot(x,y)
plt.xlim([0,30])
plt.show()
