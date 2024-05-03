import matplotlib.pyplot as plt
from scipy.io import loadmat
import numpy as np
it=10
path="Inc3HR1/"
#file="step"+str(it)+".mat"
file=(path+'step{:d}.mat'.format(it))
data1=loadmat(file)
x=data1['x'].flatten()
y=data1['y'].flatten()
phi=data1['phi']
eta_phi=data1['eta_phi']
nx,ny=np.shape(phi)
#plt.pcolor(x,y,phi.T)
plt.subplot(1, 2, 1)
#plt.contour(x, y, phi.T, cmap='viridis')
plt.pcolor(x, y, phi.T,shading='auto', cmap='viridis')
plt.axis("scaled")
plt.colorbar()
#plt.title('step {}:({:5.2f}δt)'.format(it,time/tsc))
plt.xlabel('X')
plt.ylabel('Y')
plt.subplot(1, 2, 2)
plt.plot(phi[int(nx / 2), :], y, 'r')
plt.title('Porosity Profile')
#plt.savefig(f'step{it}.jpg')
#plt.close()


