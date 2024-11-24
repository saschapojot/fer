import  numpy as np
import matplotlib.pyplot as plt


#Temperature influence lies in coefficient of second-order terms

Tc=0.9

T=1
a=-1/2*(T-Tc)


c=0.4/np.sqrt(2)/3

b=1/4

p=np.linspace(-2.5,2.5,1000)

f=a*p**2+b*p**4

np.place(f,f>0.5,np.nan)

plt.figure(figsize=(3,3))

plt.plot(p,f,'-r',linewidth=2)

plt.tick_params(axis="x",direction="in")
plt.tick_params(axis="y",direction="in")

plt.xticks([])
plt.yticks([])

plt.axhline([0],linestyle="--",color="k")
plt.axvline([0],linestyle="--",color="k")
plt.savefig("func1.png")