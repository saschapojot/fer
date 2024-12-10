import numpy as np

import matplotlib.pyplot as plt
from trial2_F_subcaseI import Hessian,dF_dz,dF_d_eta

def F(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z,omega):
    val=2*(T-Tc)*alpha2*epsyy**2 +4*theta2*epsyy**4 \
        +(2*(T-Tc)*alpha2*epsyy+8*theta2*epsyy**3)*eta \
        +(8*theta2*epsyy**2+(T-Tc)*alpha2)*eta**2 \
        +((T-Tc)*alpha1+2*beta*epsyy+gamma*epsyy**2)*z**2 \
        + omega*z*eta+(beta+gamma*epsyy)*z**2*eta \
        +4*theta2*epsyy*eta**3+theta1*z**4+theta2*eta**4

    return val


T=0.5
Tc=0.51

beta=0.7

gamma=np.sqrt(3)/2
# epsyy=-0.1

alpha1=1/3
alpha2=np.sqrt(5)

z_const=1
epsyy_const=-0.1
epsxx_const=0.2
eta_const=epsxx_const-epsyy_const
print(f"eta_const={eta_const}")
omega=-2*(beta+gamma*epsyy_const)*z_const
theta1=-((T-Tc)*alpha1+2*beta*epsyy_const+gamma*epsyy_const**2)/(2*z_const**2)
theta2=-(2*(T-Tc)*alpha2*epsyy_const+2*(T-Tc)*alpha2*eta_const+omega*z_const+(beta+gamma*epsyy_const)*z_const**2)/(8*epsyy_const**3+16*epsyy_const**2*eta_const+12*epsyy_const*eta_const**2+4*eta_const**3)

def F_2_plt(eta,z):
    # omega=-2*(beta+gamma*epsyy)*z
    # theta1=-((T-Tc)*alpha1+2*beta*epsyy*gamma*epsyy**2)/(2*z**2)
    # theta2=-(2*(T-Tc)*alpha2*epsyy+2*(T-Tc)*alpha2*eta+omega*z+(beta+gamma*epsyy)*z**2)/(8*epsyy**3+16*epsyy**2*eta+12*epsyy*eta**2+4*eta**3)

    val=F(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy_const,eta,z,omega)

    return val

H=Hessian(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy_const,eta_const,z_const,omega)
print(H)
print(dF_d_eta(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy_const,eta_const,z_const,omega))
print(dF_dz(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy_const,eta_const,z_const,omega))
eta_vec=np.linspace(-1,1,100)
z_vec=np.linspace(-2,2,100)

Eta,Z=np.meshgrid(eta_vec,z_vec)


F_PLT=F_2_plt(Eta,Z)

plt.contourf(Eta, Z, F_PLT, levels=20, cmap='coolwarm')
plt.colorbar(label='F')
plt.title('F')
plt.xlabel(r'$\eta$')
plt.ylabel('$M$')
plt.show()
plt.savefig("Ftmp.png")