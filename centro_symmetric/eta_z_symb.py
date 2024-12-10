from sympy import *
import numpy as np
from sympy.abc import alpha

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



