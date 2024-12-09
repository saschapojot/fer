import numpy as np
#subcase 1 for solving eta and z in  F
#variables

T=0.5
Tc=1
alpha1=1/3
alpha2=np.sqrt(3)

gamma=np.pi

theta1=3/10
theta2=3/5

epsyy=0.01
beta=-gamma*epsyy

coef_dF_d_eta=[4*theta2,-12*theta2*epsyy,2*((T-Tc)*alpha2+8*theta2*epsyy**2),-2*(T-Tc)*alpha2*epsyy-8*theta2*epsyy**3]

dF_d_eta_roots=np.roots(coef_dF_d_eta)


z2=(gamma*epsyy**2-(T-Tc)*alpha1)/(2*theta1)

def dF_d_eta(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z):
   val= -2*(T-Tc)*alpha2*epsyy-8*theta2*epsyy**3 \
    +2*((T-Tc)*alpha2+8*theta2*epsyy**2)*eta \
    -12*theta2*epsyy*eta**2+4*theta2*eta**3

   return val
def dF_dz(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z):
    val=2*((T-Tc)*alpha1-gamma*epsyy**2)*z+2*(beta+gamma*epsyy)*z*eta+4*theta1*z**3
    return val
def d2F_d_eta2(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z):

    val=2*((T-Tc)*alpha2+8*theta2*epsyy**2)-24*theta2*epsyy*eta+12*theta2*eta**2
    return val

def d2F_dz2(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z):
    val=2*((T-Tc)*alpha1-gamma*epsyy**2)+2*(beta+gamma*epsyy)*eta+12*theta1*z**2
    return val

def d2F_dz_d_eta(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z):
    return 2*(beta+gamma*epsyy)*z



def Hessian(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z):
    H11=d2F_d_eta2(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z)

    H12=d2F_dz_d_eta(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z)

    H21=d2F_dz_d_eta(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z)

    H22=d2F_dz2(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z)

    H=np.array([[H11,H12],[H21,H22]])
    return H


eta=dF_d_eta_roots[1]


z=np.sqrt(z2)
print(Hessian(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z))
print(dF_dz(T,Tc,alpha1,alpha2,beta,gamma,theta1,theta2,epsyy,eta,z))