import numpy as np, numba as nb
H=3600.0
num=lambda x:x/H
# ---------- parameters (from NeuronModel.nb) ----------
gl,El,gNa,ENa,SNa,V0Na,gK,EK,SK,V0K,tK=0.1,-60,3,50,0.25,-25,4,-90,0.25,-25,2
gGLU,EGLU,tGLU,gGABA,EGABA,tGABA,Ssyn,Vsp=0.1,50,30,0.1,-65,30,1,-20
gOX,EOX,tOX,tinc,tdec,Mmax=0.2,50,300,27000,3312,1
gMC,EMC,tacoa,telia,tacoA,teliA,tactMC,tinMC=0.1,-38,400,300,300,300,734.4,2030.4
basal,m=0.5,1
KK=lambda kapp,mx: kapp*basal/mx
Lx=(1-basal)/basal
kk=lambda kapp,mx: KK(kapp,mx)*Lx*mx/(1-mx)
TA,ZA,Ta,Za=KK(3.6e-9,0.1),kk(3.6e-9,0.1),KK(0.16e-9,0.9),kk(0.16e-9,0.9)
k12,k2x,sig,khx,beta,beta2,gam=num(2.30),num(0.71),0.20,num(2.354),num(1591),num(1238),1.29
th,de,s,kix,rho=num(528.63),0.38,4,num(1.03),num(90)
eps,FM,BV=num(6231.39),60000,5000
R,eta,l,EC50LEP=num(0.11),24,1,661
gINSRA,tINa,tINi,EC50INSR=0.1,300,300,661
gGHSR,tGa,tGi,EC50GHSR=0.1,300,300,500
gLEPR,tLa,tLi,EC50LEPR=0.045,300,300,16
wc,Ac=0.000073,0.045
P=np.array([gl,El,gNa,ENa,SNa,V0Na,gK,EK,SK,V0K,tK,gGLU,EGLU,tGLU,gGABA,EGABA,tGABA,Ssyn,Vsp,
  gOX,EOX,tOX,tinc,tdec,Mmax,gMC,EMC,tacoa,telia,tacoA,teliA,tactMC,tinMC,m,TA,ZA,Ta,Za,
  k12,k2x,sig,khx,beta,beta2,gam,th,de,s,kix,rho,eps,FM,BV,R,eta,l,EC50LEP,
  gINSRA,tINa,tINi,EC50INSR,gGHSR,tGa,tGi,EC50GHSR,gLEPR,tLa,tLi,EC50LEPR,wc,Ac],dtype=np.float64)

def meals(diet):
    gen=lambda h:h*H
    out=[]
    for d in (0,24,48):
        if diet=='TRE':  out+= [(gen(8+d),gen(8+d)+2340,num(2.1)),(gen(13+d),gen(13+d)+1676,num(2.1)),(gen(17+d),gen(17+d)+1400,num(2.1))]
        if diet=='SM':   out+= [(gen(8+d),gen(8+d)+2340,num(1.1)),(gen(13+d),gen(13+d)+1676,num(1.1)),(gen(17+d),gen(17+d)+1400,num(1.1))]
        if diet=='CD':   out+= [(gen(8+d),gen(8+d)+2340,num(2.1)),(gen(13+d),gen(13+d)+1676,num(2.1)),(gen(16+d),gen(16+d)+1400,num(0.8)),
                               (gen(19+d),gen(19+d)+1400,num(2.1)),(gen(20.5+d),gen(20.5+d)+1400,num(0.8))]
        if diet=='IF' and d==0: out+= [(gen(13),gen(13)+1676,num(2.1)),(gen(17),gen(17)+1400,num(2.1))]
    return np.array(out,dtype=np.float64)

@nb.njit(cache=True)
def sigm(V,S,V0): return 1.0/(1.0+np.exp(-S*(V-V0)))

@nb.njit(cache=True)
def rhs(t,y,dy,p,ml,flags):
    (gl,El,gNa,ENa,SNa,V0Na,gK,EK,SK,V0K,tK,gGLU,EGLU,tGLU,gGABA,EGABA,tGABA,Ssyn,Vsp,
     gOX,EOX,tOX,tinc,tdec,Mmax,gMC,EMC,tacoa,telia,tacoA,teliA,tactMC,tinMC,m,TA,ZA,Ta,Za,
     k12,k2x,sig,khx,beta,beta2,gam,th,de,s,kix,rho,eps,FM,BV,R,eta,l,EC50LEP,
     gINSRA,tINa,tINi,EC50INSR,gGHSR,tGa,tGi,EC50GHSR,gLEPR,tLa,tLi,EC50LEPR,wc,Ac)=p
    V,A,Pm,L,aKV,aKA,aKP,aKL,aGLUP,aGLUV,GABAA,GABAL,aOx,Mox,aMSH,AgRP,AMC4R,GL1,GL2,GHR,INS,LEP,AINSR,AGHSR,ALEPR=y
    lep_sign,ghr_paper,lepr_fix=flags[0],flags[1],flags[2]
    def Ipc(v,ak): return gl*(v-El)+gNa*sigm(v,SNa,V0Na)*(v-ENa)+gK*ak*(v-EK)
    IMC_V=gMC*AMC4R*(V-EMC); IMC_L=gMC*AMC4R*(L-EMC)
    circ=Ac*(0.97*np.sin(wc*t-np.pi/3)+0.22*np.sin(2*wc*t-np.pi/3)+0.07*np.sin(3*wc*t-np.pi/3)+0.03*np.sin(4*wc*t-np.pi/3)+0.001*np.sin(5*wc*t-np.pi/3)+1)
    dy[0]=-Ipc(V,aKV)-gGLU*aGLUP*(V-EGLU)-gGABA*GABAA*(V-EGABA)-gGABA*GABAL*(V-EGABA)-IMC_V
    dy[1]=-Ipc(A,aKA)-gGLU*aGLUV*(A-EGLU)-gOX*aOx*(A-EOX)-gGHSR*AGHSR*(A+10)-gLEPR*ALEPR*(A+90)-gINSRA*AINSR*(A+90)
    vlep = Pm if lepr_fix==1 else A
    dy[2]=-Ipc(Pm,aKP)-gGABA*GABAA*(Pm-EGABA)+gOX*aOx*(Pm-EOX)-gINSRA*AINSR*(Pm+10)-gLEPR*ALEPR*(vlep+20)
    dy[3]=-Ipc(L,aKL)-IMC_L+circ
    dy[4]=(sigm(V,SK,V0K)-aKV)/tK; dy[5]=(sigm(A,SK,V0K)-aKA)/tK
    dy[6]=(sigm(Pm,SK,V0K)-aKP)/tK; dy[7]=(sigm(L,SK,V0K)-aKL)/tK
    dy[8]=(sigm(Pm,Ssyn,Vsp)-aGLUP)/tGLU; dy[9]=(sigm(V,Ssyn,Vsp)-aGLUV)/tGLU
    dy[10]=(sigm(A,Ssyn,Vsp)-GABAA)/tGABA; dy[11]=(sigm(L,Ssyn,Vsp)-GABAL)/tGABA
    oinf=sigm(L,Ssyn,Vsp)
    dy[12]=(oinf*Mox-aOx)/tOX; dy[13]=(Mmax-Mox)/tinc-oinf*Mox/tdec
    dy[14]=1e-9*sigm(Pm,Ssyn,Vsp)/tacoa-aMSH/telia
    dy[15]=1e-8*sigm(A,Ssyn,Vsp)/tacoA-AgRP/teliA
    Ra=(aMSH+Ta*(1+AgRP/TA))/(Ta*((1+m)+((1/TA)+(m/ZA))*AgRP)+(1+(m*Ta/Za))*aMSH)
    dy[16]=Ra/tactMC-AMC4R/tinMC
    F=0.0
    for i in range(ml.shape[0]):
        if ml[i,0]<t<ml[i,1]: F=ml[i,2]
    dy[17]=-(GL1*k12)/(1+sig*GL2)+F
    dy[18]=k12/(1+sig*GL2)*GL1-k2x*GL2
    hd=(t/3600.0)%24.0
    if ghr_paper==1: b= beta if hd<6.0 else beta2
    else:            b= beta2 if (3.0<=hd<6.0 and t>=3*3600.0) else beta
    dy[19]=b/(1+gam*GL2)-khx*GHR
    dy[20]=rho+th*GL2**s/(de**s+GL2**s)-kix*INS
    dy[21]=(FM*(R*eta)/(1+(INS/EC50LEP)**(lep_sign*l))-eps*LEP)/BV
    dy[22]=(INS/(EC50INSR+INS))/tINa-AINSR/tINi
    dy[23]=(GHR/(EC50GHSR+GHR))/tGa-AGHSR/tGi
    dy[24]=(LEP/(EC50LEPR+LEP))/tLa-ALEPR/tLi

@nb.njit(cache=True)
def run(p,ml,flags,T,dt,rec_every):
    y=np.array([-40.,-60,-60,-60,0,0,0,0,0,0,0,0,0,1,0,0,1.38,0,0.01,670,79,22,0.11,0.6,0.6])
    n=y.size; k1=np.empty(n);k2=np.empty(n);k3=np.empty(n);k4=np.empty(n);yt=np.empty(n)
    N=int(T/dt); nrec=N//rec_every+1
    rec=np.empty((nrec,n)); spikes=np.zeros((4,200000)); ns=np.zeros(4,np.int64)
    prev=y[:4].copy(); t=0.0; r=0
    for i in range(N):
        if i%rec_every==0: rec[r]=y; r+=1
        rhs(t,y,k1,p,ml,flags)
        for j in range(n): yt[j]=y[j]+0.5*dt*k1[j]
        rhs(t+0.5*dt,yt,k2,p,ml,flags)
        for j in range(n): yt[j]=y[j]+0.5*dt*k2[j]
        rhs(t+0.5*dt,yt,k3,p,ml,flags)
        for j in range(n): yt[j]=y[j]+dt*k3[j]
        rhs(t+dt,yt,k4,p,ml,flags)
        for j in range(n): y[j]+=dt*(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6
        t+=dt
        for q in range(4):
            if prev[q]<-20.0 and y[q]>=-20.0 and ns[q]<200000:
                spikes[q,ns[q]]=t; ns[q]+=1
            prev[q]=y[q]
    return rec[:r],spikes,ns
