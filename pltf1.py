import  numpy as np
import matplotlib.pyplot as plt



#C3-breaking ferroelectricity
t=np.linspace(0,2*np.pi,1000) #angle theta

r=np.linspace(0,2.3,1000)

T,R=np.meshgrid(t,r)

f=0.1*(19-20)*R**2/2+1*R**4/4-0.4*R**3*np.cos(3*T)/3
np.place(f,f>0.0,np.nan)
fig=plt.figure()
ax=plt.subplot(1,1,1,projection="polar")
contour=ax.contour(t,r,f,cmap="coolwarm",levels=20)
# plt.contour(t,r,f)
cbar=plt.colorbar(contour,ax=ax)
plt.savefig("func2.png")



