# Modelo hipotalámico en HORAS, parámetros convertidos explícitamente (espejo de HypothalamusModel.wl)
import numpy as np, numba as nb
S=1/3600.0   # 1 s en horas
Cm=S         # capacitancia efectiva: C dV/dt con t en h (C=1 s en el original)
p=dict(gl=0.1,El=-60,gNa=3,ENa=50,SNa=0.25,V0Na=-25,gK=4,EK=-90,SK=0.25,V0K=-25,tK=2*S,
 gGLU=0.1,EGLU=50,tGLU=30*S,gGABA=0.1,EGABA=-65,tGABA=30*S,Ssyn=1,Vsp=-20,
 gOX=0.2,EOX=50,tOX=300*S,tinc=27000*S,tdec=3312*S,Mmax=1,
 gMC=0.1,EMC=-38,tacoa=400*S,telia=300*S,tacoA=300*S,teliA=300*S,tactMC=734.4*S,tinMC=2030.4*S,m=1,
 k12=2.30,k2x=0.71,sig=0.20,khx=2.354,beta=1591,beta2=1238,gam=1.29,
 th=528.63,de=0.38,s=4,kix=1.03,rho=90,eps=6231.39,FM=60000,BV=5000,R=0.11,eta=24,l=1,EC50LEP=661,
 gINSRA=0.1,tINa=300*S,tINi=300*S,EC50INSR=661,gGHSR=0.1,tGa=300*S,tGi=300*S,EC50GHSR=500,
 gLEPR=0.045,tLa=300*S,tLi=300*S,EC50LEPR=16,wc=0.000073*3600,Ac=0.045)
basal=0.5; KK=lambda k,mx:k*basal/mx; kk=lambda k,mx:KK(k,mx)*((1-basal)/basal)*mx/(1-mx)
p.update(TA=KK(3.6e-9,0.1),ZA=kk(3.6e-9,0.1),Ta=KK(0.16e-9,0.9),Za=kk(0.16e-9,0.9))
KEYS=list(p.keys()); P=np.array([p[k] for k in KEYS])
DIETS={'TRE':[(8,0.65,2.1),(13,1676/3600,2.1),(17,1400/3600,2.1)],
       'SM': [(8,0.65,1.1),(13,1676/3600,1.1),(17,1400/3600,1.1)],
       'CD': [(8,0.65,2.1),(13,1676/3600,2.1),(16,1400/3600,0.8),(19,1400/3600,2.1),(20.5,1400/3600,0.8)],
       'IF': [(13,1676/3600,2.1),(17,1400/3600,2.1)]}
@nb.njit(cache=True)
def sg(V,S_,V0): return 1/(1+np.exp(-S_*(V-V0)))
@nb.njit(cache=True)
def rhs(t,y,dy,q,ml):
    (gl,El,gNa,ENa,SNa,V0Na,gK,EK,SK,V0K,tK,gGLU,EGLU,tGLU,gGABA,EGABA,tGABA,Ssyn,Vsp,gOX,EOX,tOX,tinc,tdec,Mmax,
     gMC,EMC,tacoa,telia,tacoA,teliA,tactMC,tinMC,m,k12,k2x,sig,khx,beta,beta2,gam,th,de,s,kix,rho,eps,FM,BV,R,eta,l,EC50LEP,
     gINSRA,tINa,tINi,EC50INSR,gGHSR,tGa,tGi,EC50GHSR,gLEPR,tLa,tLi,EC50LEPR,wc,Ac,TA,ZA,Ta,Za)=q
    Cm=1/3600.0
    V,A,Pm,L,aKV,aKA,aKP,aKL,aGLUP,aGLUV,GABAA,GABAL,aOx,Mox,aMSH,AgRP,AMC4R,GL1,GL2,GHR,INS,LEP,AINSR,AGHSR,ALEPR=y
    def Ipc(v,ak): return gl*(v-El)+gNa*sg(v,SNa,V0Na)*(v-ENa)+gK*ak*(v-EK)
    circ=Ac*(0.97*np.sin(wc*t-np.pi/3)+0.22*np.sin(2*wc*t-np.pi/3)+0.07*np.sin(3*wc*t-np.pi/3)+0.03*np.sin(4*wc*t-np.pi/3)+0.001*np.sin(5*wc*t-np.pi/3)+1)
    dy[0]=(-Ipc(V,aKV)-gGLU*aGLUP*(V-EGLU)-gGABA*GABAA*(V-EGABA)-gGABA*GABAL*(V-EGABA)-gMC*AMC4R*(V-EMC))/Cm
    dy[1]=(-Ipc(A,aKA)-gGLU*aGLUV*(A-EGLU)-gOX*aOx*(A-EOX)-gGHSR*AGHSR*(A+10)-gLEPR*ALEPR*(A+90)-gINSRA*AINSR*(A+90))/Cm
    dy[2]=(-Ipc(Pm,aKP)-gGABA*GABAA*(Pm-EGABA)+gOX*aOx*(Pm-EOX)-gINSRA*AINSR*(Pm+10)-gLEPR*ALEPR*(Pm+20))/Cm
    dy[3]=(-Ipc(L,aKL)-gMC*AMC4R*(L-EMC)+circ)/Cm
    dy[4]=(sg(V,SK,V0K)-aKV)/tK; dy[5]=(sg(A,SK,V0K)-aKA)/tK; dy[6]=(sg(Pm,SK,V0K)-aKP)/tK; dy[7]=(sg(L,SK,V0K)-aKL)/tK
    dy[8]=(sg(Pm,Ssyn,Vsp)-aGLUP)/tGLU; dy[9]=(sg(V,Ssyn,Vsp)-aGLUV)/tGLU
    dy[10]=(sg(A,Ssyn,Vsp)-GABAA)/tGABA; dy[11]=(sg(L,Ssyn,Vsp)-GABAL)/tGABA
    oi=sg(L,Ssyn,Vsp); dy[12]=(oi*Mox-aOx)/tOX; dy[13]=(Mmax-Mox)/tinc-oi*Mox/tdec
    dy[14]=1e-9*sg(Pm,Ssyn,Vsp)/tacoa-aMSH/telia; dy[15]=1e-8*sg(A,Ssyn,Vsp)/tacoA-AgRP/teliA
    Ra=(aMSH+Ta*(1+AgRP/TA))/(Ta*((1+m)+((1/TA)+(m/ZA))*AgRP)+(1+(m*Ta/Za))*aMSH)
    dy[16]=Ra/tactMC-AMC4R/tinMC
    hd=t%24.0; F=0.0
    for i in range(ml.shape[0]):
        if ml[i,0]<hd<ml[i,0]+ml[i,1]: F=ml[i,2]
    dy[17]=-(GL1*k12)/(1+sig*GL2)+F; dy[18]=k12/(1+sig*GL2)*GL1-k2x*GL2
    b=beta2 if 3.0<=hd<6.0 else beta
    dy[19]=b/(1+gam*GL2)-khx*GHR; dy[20]=rho+th*GL2**s/(de**s+GL2**s)-kix*INS
    dy[21]=(FM*R*eta/(1+(INS/EC50LEP)**l)-eps*LEP)/BV
    dy[22]=(INS/(EC50INSR+INS))/tINa-AINSR/tINi; dy[23]=(GHR/(EC50GHSR+GHR))/tGa-AGHSR/tGi; dy[24]=(LEP/(EC50LEPR+LEP))/tLa-ALEPR/tLi
@nb.njit(cache=True)
def run(q,ml,T,dt,rec_every):
    y=np.array([-40.,-60,-60,-60,0,0,0,0,0,0,0,0,0,1,0,0,1.38,0,0.01,670,79,22,0.11,0.6,0.6])
    n=25;k1=np.empty(n);k2=np.empty(n);k3=np.empty(n);k4=np.empty(n);yt=np.empty(n)
    N=int(round(T/dt)); rec=np.empty((N//rec_every+1,n)); sp=np.zeros((4,400000)); ns=np.zeros(4,np.int64)
    prev=y[:4].copy(); t=0.0; r=0
    for i in range(N):
        if i%rec_every==0: rec[r]=y; r+=1
        rhs(t,y,k1,q,ml)
        for j in range(n): yt[j]=y[j]+0.5*dt*k1[j]
        rhs(t+0.5*dt,yt,k2,q,ml)
        for j in range(n): yt[j]=y[j]+0.5*dt*k2[j]
        rhs(t+0.5*dt,yt,k3,q,ml)
        for j in range(n): yt[j]=y[j]+dt*k3[j]
        rhs(t+dt,yt,k4,q,ml)
        for j in range(n): y[j]+=dt*(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6
        t+=dt
        for qq in range(4):
            if prev[qq]<-20.0 and y[qq]>=-20.0: sp[qq,ns[qq]]=t; ns[qq]+=1
            prev[qq]=y[qq]
    return rec[:r],sp,ns
def simulate(diet,hours,dt=0.01/3600):
    rec,sp,ns=run(P,np.array(DIETS[diet],dtype=np.float64),float(hours),dt,1000)
    return rec,[sp[i,:ns[i]].copy() for i in range(4)]
