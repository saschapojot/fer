import numpy as np

#trial 2 of F, subcase I

T=0.5
Tc=0.51

beta=0.7

gamma=np.sqrt(3)/2
epsyy=-0.1
epsxx=0.2
alpha1=1/3
alpha2=np.sqrt(5)
eta=epsxx-epsyy
z=-1
omega=-2*(beta+gamma*epsyy)*z

theta1=-((T-Tc)*alpha1+2*beta*epsyy*gamma*epsyy**2)/(2*z**2)

theta2=-(2*(T-Tc)*alpha2*epsyy+2*(T-Tc)*alpha2*eta+omega*z+(beta+gamma*epsyy)*z**2)/(8*epsyy**3+16*epsyy**2*eta+12*epsyy*eta**2+4*eta**3)


def d2F_d_eta2(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z,omega):
    val=2*(8*theta2*epsyy**2+(T-Tc)*alpha2)+24*theta2*epsyy*eta+12*theta2*eta**2
    return val

def d2F_dz2(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z,omega):
    val=2*((T-Tc)*alpha1+2*beta*epsyy+gamma*epsyy**2)+2*(beta+gamma*epsyy)*eta+12*theta1*z**2
    return val


def d2F_dz_d_eta(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z,omega):
    val=omega+2*(beta+gamma*epsyy)*z

    return val


def Hessian(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z,omega):
    H11=d2F_d_eta2(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z,omega)

    H12=d2F_dz_d_eta(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z,omega)

    H21=d2F_dz_d_eta(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z,omega)

    H22=d2F_dz2(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z,omega)

    H=np.array([[H11,H12],[H21,H22]])

    return H

H=Hessian(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z,omega)
print(H)