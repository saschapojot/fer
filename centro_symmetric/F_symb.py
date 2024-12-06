from sympy import *


T,Tc,alpha1,alpha2=symbols("T,Tc,alpha1,alpha2",cls=Symbol,real=True)

beta,gamma,theta1,theta2=symbols("beta,gamma,theta1,theta2",cls=Symbol,real=True)

z,eta=symbols(r"z,eta",cls=Symbol,real=True)

epsyy=symbols("epsilon_yy",cls=Symbol,real=True)

epsxx=eta+epsyy


F_original=(T-Tc)*(alpha1*z**2+alpha2*(epsxx**2+epsyy**2))\
         +beta*z*(epsxx-epsyy)+gamma*z*epsxx*epsyy\
        + theta1*z**4+theta2*(epsxx**2+epsyy**2)**2

F=2*(T-Tc)*alpha2*epsyy**2+4*theta2*epsyy**4\
    +(2*(T-Tc)*alpha2*epsyy+8*theta2*epsyy**3)*eta\
    +gamma*epsyy**2*z\
    +(T-Tc)*alpha1*z**2+(beta+gamma*epsyy)*z*eta+((T-Tc)*alpha2+8*theta2*epsyy**2)*eta**2\
    +4*theta2*epsyy*eta**3+theta1*z**4+theta2*eta**4




dF_d_eta=2*(T-Tc)*alpha2*epsyy+ 8*theta2*epsyy**3+(beta+gamma*epsyy)*z\
    +2*((T-Tc)*alpha2+8*theta2*epsyy**2)*eta+12*theta2*epsyy*eta**2+4*theta2*eta**3

dF_dz=gamma*epsyy**2+2*(T-Tc)*alpha1*z+(beta+gamma*epsyy)*eta+4*theta1*z**3

d2F_dz_d_eta=beta+gamma*epsyy

d2F_d_eta2=2*((T-Tc)*alpha2+8*theta2*epsyy**2)+24*theta2*epsyy*eta+12*theta2*eta**2

d2F_dz2=2*(T-Tc)*alpha1+12*theta1*z**2

