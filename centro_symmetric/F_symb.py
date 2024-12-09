from sympy import *
import numpy as np

T,Tc,alpha1,alpha2=symbols("T,Tc,alpha1,alpha2",cls=Symbol,real=True)

beta,gamma,theta1,theta2=symbols("beta,gamma,theta1,theta2",cls=Symbol,real=True)

z,eta=symbols(r"z,eta",cls=Symbol,real=True)

epsyy=symbols("epsilon_yy",cls=Symbol,real=True)
T=0.5
Tc=1
alpha1=1/3
alpha2=np.sqrt(3)

gamma=np.pi

theta1=3/10
theta2=3/5

epsyy=0.01
beta=-gamma*epsyy
eta=1
z2=(gamma*epsyy**2-(T-Tc)*alpha1)/(2*theta1)
z=np.sqrt(z2)
epsxx=eta-epsyy


F_original=(T-Tc)*(alpha1*z**2+alpha2*(epsxx**2+epsyy**2))\
         +beta*z**2*(epsxx+epsyy)+gamma*z**2*epsxx*epsyy\
        + theta1*z**4+theta2*(epsxx**2+epsyy**2)**2


F=4*theta2*epsyy**4+2*(T-Tc)*alpha2*epsyy**2 + (-2*(T-Tc)*alpha2*epsyy-8*theta2*epsyy**3)*eta \
    +((T-Tc)*alpha2+8*theta2*epsyy**2)*eta**2\
    +((T-Tc)*alpha1-gamma*epsyy**2)*z**2\
    +(beta+gamma*epsyy)*z**2*eta\
    -4*theta2*epsyy*eta**3+theta1*z**4+theta2*eta**4


dF_d_eta=-2*(T-Tc)*alpha2*epsyy-8*theta2*epsyy**3\
        +2*((T-Tc)*alpha2+8*theta2*epsyy**2)*eta\
        +(beta+gamma*epsyy)*z**2-12*theta2*epsyy*eta**2+4*theta2*eta**3

dF_dz=2*((T-Tc)*alpha1-gamma*epsyy**2)*z+2*(beta+gamma*epsyy)*z*eta+4*theta1*z**3
d2F_dzd_eta=2*(beta+gamma*epsyy)*z

d2F_d_eta2=2*((T-Tc)*alpha2+8*theta2*epsyy**2)-24*theta2*epsyy*eta+12*theta2*eta**2

d2F_dz2=2*((T-Tc)*alpha1-gamma*epsyy**2)+2*(beta+gamma*epsyy)*eta+12*theta1*z**2


pprint(d2F_dz2)








