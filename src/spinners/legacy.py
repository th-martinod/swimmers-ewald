"""Faithful numpy port of base_sums.rkt (production version), kept for regression only.

Reproduces the published behaviour, including its known defects:
- Fourier part evaluated with Q rotated by +30 deg for 'hex' (randdirf) and projected on b1,
- lattice 'A' factor in the velocity sums, spurious h^3/pi^3 prefactor,
- real-space prefactor 4/pi too large (hidden in production because alpha = 1 << h^2).
Do NOT use for new results. See docs/LEGACY_DEFECTS.md.
"""
import numpy as np
from scipy.special import erfc
pi=np.pi
def latt_k(latt,h):
    if latt=="scl": B=np.array([[1,0,0],[0,1,0],[0,0,1.]])/h
    if latt=="bcl": B=np.array([[1,1,0],[0,1,1],[1,0,1.]])/h
    if latt=="fcl": B=np.array([[1,1,-1],[-1,1,1],[1,-1,1.]])/h
    if latt=="hex": B=np.array([[1/h,1/(np.sqrt(3)*h),0],[0,2/(np.sqrt(3)*h),0],[0,0,1/h]])
    return B
def latt_r(latt,h):
    if latt=="scl": A=np.eye(3)*h
    if latt=="bcl": A=np.array([[1,1,-1],[-1,1,1],[1,-1,1.]])*h/2
    if latt=="fcl": A=np.array([[1,1,0],[0,1,1],[1,0,1.]])*h/2
    if latt=="hex": A=np.array([[h,0,0],[-h/2,np.sqrt(3)*h/2,0],[0,0,h]])
    return A
def V0f(latt,h): return {"scl":h**3,"bcl":h**3/2,"fcl":h**3/4,"hex":np.sqrt(3)*h**3/2}[latt]
def idx(n,n3):
    r=np.arange(-n,n+1); r3=np.arange(-n3,n3+1)
    I,J,K=np.meshgrid(r,r,r3,indexing='ij'); I,J,K=I.ravel(),J.ravel(),K.ravel()
    m=~((I==0)&(J==0)&(K==0)); return I[m],J[m],K[m]
def Bdot(k,x):  # matB(k).x for rows of k, x shape (N,3) or (3,)
    k2=(k*k).sum(1)
    return k*((k*x).sum(1)/k2**2)[:,None] - x/k2[:,None]
def phi1d2(r): return np.exp(-r)/r + np.sqrt(pi)/(2*r**1.5)*erfc(np.sqrt(r))
def Adot(alpha,r,x):  # matA(alpha,r).x
    r2=(r*r).sum(1)
    t1=2*alpha**-0.5*np.exp(-pi*r2/alpha)
    t2=pi*alpha**-1.5*phi1d2(pi*r2/alpha)
    return (t1-t2*r2)[:,None]*x - t2[:,None]*r*((r*x).sum(1))[:,None]
def conf_fourier(q,nv,b,h,alpha,n,latt,dim):
    n3=0 if dim=="2D" else n
    I,J,K=idx(n,n3); k=np.outer(I,latt_k(latt,h)[0])+np.outer(J,latt_k(latt,h)[1])+np.outer(K,latt_k(latt,h)[2])
    nq=np.linalg.norm(q); k2=(k*k).sum(1); qxn=np.cross(q,nv)
    w=(1+pi*alpha*k2)*np.exp(-pi*alpha*k2)*np.sin(2*pi*b*(k@qxn))**2*np.sin(pi*(k@q))**2
    s=Bdot(k,np.tile(q,(len(k),1))*w[:,None]).sum(0)
    return -1/(pi**2*V0f(latt,h))*h**3/pi**3/nq*s
def conf_real(q,nv,b,h,alpha,n,latt,dim):
    n3=0 if dim=="2D" else n
    I,J,K=idx(n,n3); A=latt_r(latt,h); r=np.outer(I,A[0])+np.outer(J,A[1])+np.outer(K,A[2])
    nq=np.linalg.norm(q); e=np.cross(q,nv)*b/nq
    offs=[(-q,-2),(q,-2),(-2*e,-2),(-q-2*e,1),(q-2*e,1),(2*e,-2),(-q+2*e,1),(q+2*e,1),(0*q,4)]
    s=sum(w*Adot(alpha,-r+o,np.tile(q,(len(r),1))).sum(0) for o,w in offs)
    return -1/(4*pi**2)*h**3/pi**3/nq*s
def rot_fourier(q,nv,b,h,alpha,n,latt,dim,useA=True):
    n3=0 if dim=="2D" else n
    I,J,K=idx(n,n3); kb=latt_k(latt,h); k=np.outer(I,kb[0])+np.outer(J,kb[1])+np.outer(K,kb[2])
    nq=np.linalg.norm(q); k2=(k*k).sum(1); qxn=np.cross(q,nv)
    Af=Afac(latt,I,J,K) if useA else 1.0
    w=Af*(1+pi*alpha*k2)*np.exp(-pi*alpha*k2)*np.sin(2*pi*b*(k@qxn))*np.sin(pi*(k@q))
    s=Bdot(k,np.tile(q,(len(k),1))*w[:,None]).sum(0)
    return -1/(pi**2*V0f(latt,h))*h**3/pi**3/nq*s
def rot_real(q,nv,b,h,alpha,n,latt,dim,useA=True):
    n3=0 if dim=="2D" else n
    I,J,K=idx(n,n3); A=latt_r(latt,h); r=np.outer(I,A[0])+np.outer(J,A[1])+np.outer(K,A[2])
    nq=np.linalg.norm(q); e=np.cross(q,nv)*b/nq
    Af=Afac(latt,I,J,K) if useA else np.ones(len(I))
    terms=[(q/2+e,1),(q/2-e,-1),(-q/2+e,-1),(-q/2-e,1)]
    s=(Af[:,None]*sum(w*Adot(alpha,-r+o,np.tile(q,(len(r),1))) for o,w in terms)).sum(0)
    return 1/(4*pi**2)*h**3/pi**3/nq*s
def Afac(latt,I,J,K):
    with np.errstate(all='ignore'):
        if latt=="scl":
            z=(I==0)|(J==0)|(K==0); v=np.sin(pi*I)*np.sin(pi*J)*np.sin(pi*K)/(I*J*K)
        if latt=="hex":
            z=(I==0)|(I+2*J==0)|(K==0); v=np.sqrt(3)*np.sin(pi*I)*np.sin(pi*(I+2*J)/np.sqrt(3))*np.sin(pi*K)/(I*(I+2*J)*K)
    return np.where(z,1.0,v)
def dirs(latt,part):
    if part=="F": return {"scl":[1,0,0],"bcl":[1,1,0],"fcl":[1,1,-1],"hex":[1,1/np.sqrt(3),0]}[latt]
    return {"scl":[1,0,0],"bcl":[1,1,-1],"fcl":[1,1,0],"hex":[1,0,0]}[latt]
def nz(v):
    v=np.array(v,float); n=np.linalg.norm(v); return v/n if n>0 else 0*v
def conf_total(th,latt,dim,h,alpha,nF,nR,b=0.5,l0=1.0,bug=True):
    nv=np.array([0,0,1.])
    q=l0*np.array([np.cos(th),np.sin(th),0])
    if latt=="hex" and bug:
        rf=np.array([(3*np.cos(th)-np.sqrt(3)*np.sin(th))/3,(np.sqrt(3)*np.cos(th)+3*np.sin(th))/3,0])
        qf=l0*nz(rf); dF=dirs("hex","F")
    else: qf=q; dF=dirs(latt,"R")
    vF=conf_fourier(qf,nv,b,h,alpha,nF,latt,dim); vR=conf_real(q,nv,b,h,alpha,nR,latt,dim)
    aF=np.cross(qf,vF)@nz(np.cross(dF,qf))/np.linalg.norm(qf)
    aR=np.cross(q,vR)@nz(np.cross(dirs(latt,"R"),q))/np.linalg.norm(q)
    return aF+aR, aF, aR
