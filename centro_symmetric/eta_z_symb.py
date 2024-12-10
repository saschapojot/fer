from sympy import *
import numpy as np
from sympy.abc import alpha

#trial 2, diagonal mirror symmetry
T,Tc,alpha1,alpha2=symbols("T,Tc,alpha1,alpha2",cls=Symbol,real=True)

beta,gamma,theta1,theta2=symbols("beta,gamma,theta1,theta2",cls=Symbol,real=True)

z,eta=symbols(r"z,eta",cls=Symbol,real=True)

epsyy=symbols("epsilon_yy",cls=Symbol,real=True)

omega=symbols("omega",cls=Symbol,real=True)
# epsxx=symbols("epsilon_xx",cls=Symbol,real=True)
epsxx=eta+epsyy

F_original=(T-Tc)*alpha1*z**2 + (T-Tc)*alpha2*(epsxx**2+epsyy**2)\
    +omega*z*(epsxx-epsyy) +beta*z**2*(epsxx+epsyy)+ gamma*z**2*epsxx*epsyy\
           +theta1*z**4+theta2*(epsxx**2+epsyy**2)**2



F=2*(T-Tc)*alpha2*epsyy**2 +4*theta2*epsyy**4\
    +(2*(T-Tc)*alpha2*epsyy+8*theta2*epsyy**3)*eta\
    +(8*theta2*epsyy**2+(T-Tc)*alpha2)*eta**2\
    +((T-Tc)*alpha1+2*beta*epsyy+gamma*epsyy**2)*z**2\
    + omega*z*eta+(beta+gamma*epsyy)*z**2*eta\
    +4*theta2*epsyy*eta**3+theta1*z**4+theta2*eta**4



dF_d_eta=2*(T-Tc)*alpha2*epsyy+8*theta2*epsyy**3\
    +2*(8*theta2*epsyy**2+(T-Tc)*alpha2)*eta\
    +omega*z+(beta+gamma*epsyy)*z**2\
    +12*theta2*epsyy*eta**2+4*theta2*eta**3

dF_dz=2*((T-Tc)*alpha1+2*beta*epsyy+gamma*epsyy**2)*z+omega*eta+2*(beta+gamma*epsyy)*z*eta+4*theta1*z**3

d2F_dz_d_eta=omega+2*(beta+gamma*epsyy)*z

d2F_d_eta2=2*(8*theta2*epsyy**2+(T-Tc)*alpha2)+24*theta2*epsyy*eta+12*theta2*eta**2

d2F_dz2=2*((T-Tc)*alpha1+2*beta*epsyy+gamma*epsyy**2)+2*(beta+gamma*epsyy)*eta+12*theta1*z**2

tmp=diff(F,(z,2))-d2F_dz2
pprint(expand(tmp))