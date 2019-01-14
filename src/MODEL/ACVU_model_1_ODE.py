import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate

# opening of operator
# 1) O* <--> O: fO,bO
# binding of HLH-2
# 2) O+H <--> OH: fOH, bOH
# production of HLH-2
# 3) 0 <--> H: fH, bH+bH*S^n/(K^n+S^n)
# production and degradation of Delta mRNA
# 4) OH --> OH + M: fM
# 5) M-->0: bM
# production and degradation of Delta
# 6) M --> M + D: fD
# 7) D --> 0: bD
# Notch signaling strength in cell i
# S_i = phi D_j

# model B: assume 3) in steady state 


def alpha_fun(S):
    return KH*(K**n+S**n)/(K**n+(1+b)*S**n)

def ODE_fun(y,t):

    O=y[[0,4]]
    OH=y[[1,5]]
    M=y[[2,6]]
    D=y[[3,7]]
    
    dOdt=np.zeros(2)
    dOHdt=np.zeros(2)
    dMdt=np.zeros(2)
    dDdt=np.zeros(2)
    
    S=phi*D
    
    c=[1,0]
    for i in range(0,2):
        j=c[i]
        dOdt[i]=fO*(1-O[i]-OH[i]) - bO*O[i] - fOH*O[i]*alpha_fun(S[j]) + bOH*OH[i]
        dOHdt[i]=fOH*O[i]*alpha_fun(S[j]) - bOH*OH[i]
        dMdt[i]=fM*OH[i] - bM*M[i]
        dDdt[i]=fD*M[i] - bD*D[i]
        
    dydt=np.array([dOdt[0],dOHdt[0],dMdt[0],dDdt[0], dOdt[1],dOHdt[1],dMdt[1],dDdt[1]])
    
    return (dydt)

n=2
KH=10.
KOH=1.
KO=0.1
K=4e2
b=5e3
KD=10.
KM=30.
phi=0.0

# x corresponds to [OH]
x=np.linspace(0,1,100)

S=phi*KD*KM*x
a=KOH*alpha_fun(S)
y=a*KO/(1+KO+a*KO)

plt.clf()
plt.subplot(2,1,1)
plt.plot(x,y,y,x)

fO=1e-2
bO=fO/KO
fOH=1e0
bOH=fOH/KOH
fM=1e0
bM=fM/KM
fD=1e-2
bD=fD/KD
    
y0=np.array([0.0,0.0,0.0,0.0, 0.1,0.0,0.0,0.0])

t=np.linspace(0,5e3,1e4)
y=integrate.odeint(ODE_fun,y0,t)

plt.plot(y[:,1],y[:,5])

plt.subplot(2,1,2)
plt.plot(t,y[:,2],'-k',t,y[:,3],'-b')
plt.plot(t,y[:,6],'--k',t,y[:,7],'--b')

print(KM*KOH*alpha_fun(0)*KO/(1+(1+KOH*alpha_fun(0))*KO), y[-1,6])

plt.show()
