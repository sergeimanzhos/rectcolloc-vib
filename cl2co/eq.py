import numpy as np
from scipy.optimize import minimize
from pes import NNPES, rij, REQ
P=NNPES('fitNN100NNs4.dat')
AMU2ME=1822.888486209; H2CM=219474.6313705
m=np.array([12.0,15.99491461957,34.968852682,34.968852682])*AMU2ME
def Xc(p):
    rCO,rCCl,a=p
    return np.array([[0,0,0],[0,0,rCO],[0,rCCl*np.sin(a),-rCCl*np.cos(a)],[0,-rCCl*np.sin(a),-rCCl*np.cos(a)]],float)
Vi=lambda p: P(rij(Xc(p)))[0]
a0=np.arccos((2*REQ[1]**2-REQ[5]**2)/(2*REQ[1]**2))/2
r=minimize(Vi,[REQ[0],REQ[1],a0],method='Nelder-Mead',options={'xatol':1e-10,'fatol':1e-11,'maxiter':40000})
r=minimize(Vi,r.x,method='Nelder-Mead',options={'xatol':1e-11,'fatol':1e-12,'maxiter':40000})
X0=Xc(r.x); X0=X0-(m[:,None]*X0).sum(0)/m.sum(); Vmin=r.fun
V=lambda x: P(rij(x.reshape(-1,4,3)))
h=2e-3; x0=X0.ravel(); n=12; H=np.zeros((n,n))
E=lambda dx: V(x0+dx)[0]
for i in range(n):
    for j in range(i,n):
        if i==j:
            e=[E(k*h*np.eye(n)[i]) for k in (-2,-1,0,1,2)]
            H[i,i]=(-e[0]+16*e[1]-30*e[2]+16*e[3]-e[4])/(12*h*h)
        else:
            ei,ej=np.eye(n)[i]*h,np.eye(n)[j]*h
            H[i,j]=H[j,i]=(E(ei+ej)-E(ei-ej)-E(-ei+ej)+E(-ei-ej))/(4*h*h)
H/=H2CM   # hartree/bohr^2
print('Vmin %.4f cm-1  rCO %.5f A rCCl %.5f A  ClCCl %.3f deg'%(Vmin,r.x[0]/1.8897261254535,r.x[1]/1.8897261254535,np.degrees(2*r.x[2])))
np.savez('eq.npz',X0=X0,m=m,Hess=H,Vmin=Vmin)
