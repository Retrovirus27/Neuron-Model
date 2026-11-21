import pickle,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
L=lambda d,h,n: pickle.load(open(f'res_{d}_{h}_{n}.pkl','rb'))
ins=np.array([[0,120.83],[1,113.24],[2,106.86],[3,103.65],[4,90.16],[5,79.23],[6,79.44],[7,89.17],[8,89.17],[8.5,71.34],[9.5,391.79],[10,518.51],[10.5,502.54],[11,373.14],[11.5,331.01],[12,214.13],[12.5,168.69],[13,384.29],[13.5,341.9],[14.5,248.41],[15,244.37],[15.5,154.12],[16,199.26],[16.5,155.32],[17,170.29],[17.5,126.07],[18,144.55],[18.5,280.14],[19,382.17],[19.5,393.84],[20,349.66],[20.5,292.61],[21,222.61],[21.5,188.44],[22,159.02],[23,129.76]])
lep=np.array([[0,22.8286],[1,23.4473],[2,23.503],[3,23.5769],[4,23.9776],[5,23.9425],[6,23.2717],[7,22.6371],[8,22.5111],[8.5,21.0934],[9,19.4697],[9.5,17.1366],[10,16.0502],[10.5,15.6503],[11,15.9026],[11.5,16.3265],[12,16.6017],[12.5,17.3575],[13,17.0149],[13.5,17.9537],[14,16.3524],[14.5,17.7032],[15,17.7882],[15.5,18.0114],[16,18.506],[16.5,19.7699],[17,19.8795],[17.5,20.316],[18,20.7525],[18.5,21.298],[19,20.6447],[19.5,19.6463],[20,19.774],[20.5,20.0834],[21,20.2112],[22,21.1926],[23,21.393]]); lep[:,1]/=2.0043
ghr=np.array([[0,625.29],[1,650.48],[2,626.89],[3,613.66],[4,542.79],[5,504.43],[6,501.48],[7,600.55],[8,704.07],[8.5,525.29],[9,321.37],[9.5,343.57],[10,352.47],[10.5,345.11],[11,389.49],[11.5,450.11],[12,648.17],[12.5,689.58],[13,432.47],[13.5,451.71],[14,422.18],[14.5,408.92],[15,502.05],[15.5,484.35],[16,509.51],[16.5,540.58],[17,691.34],[17.5,704.68],[18,639.68],[18.5,514.1],[19,515.61],[19.5,450.62],[20,492.03],[20.5,475.81],[21,602.93],[22,608.91],[23,608.97]]); ghr[:,1]/=3.3709
def hor(r): 
    t=np.arange(len(r['rec']))*10/3600; R=r['rec']; return t,R[:,19]/3.3709,R[:,20],R[:,21]/2.0043
def rmse(t,y,d): return np.sqrt(np.mean((np.interp(d[:,0],t,y)-d[:,1])**2))
def rate(r,q,h): return np.histogram(r['sp'][q]/3600,bins=np.arange(0,h+0.01,0.5))[0]
COL={'publicado':'k','LEPR_POMC_corregido':'tab:red','leptina_como_texto':'tab:blue','ghrelina_como_texto':'tab:green'}
LAB={'publicado':'Publicado (código del paper)','LEPR_POMC_corregido':'Corrección I_LEPR en POMC','leptina_como_texto':'Leptina como Ec. 9 del texto','ghrelina_como_texto':'Ghrelina como texto (β1 0–6 h)'}
out=[]
# ---- Fig 1 hormones
f,ax=plt.subplots(3,1,figsize=(7,7),sharex=True)
for nm in ['publicado','leptina_como_texto','ghrelina_como_texto']:
    t,g,i,lp=hor(L('TRE',24,nm))
    for a,y,d in zip(ax,(g,i,lp),(ghr,ins,lep)): a.plot(t,y,color=COL[nm],lw=1.6,label=LAB[nm])
    out.append((nm,rmse(t,g,ghr),rmse(t,i,ins),rmse(t,lp,lep)))
for a,d,lb in zip(ax,(ghr,ins,lep),('Ghrelina (pM)','Insulina (pM)','Leptina (pM)')):
    a.plot(d[:,0],d[:,1],'o',color='gray',ms=4,label='Datos Cummings 2001'); a.set_ylabel(lb)
    for x in (8,13,17): a.axvline(x,color='0.85',zorder=0)
ax[0].legend(fontsize=7,loc='lower left'); ax[-1].set_xlabel('Tiempo (h)'); ax[-1].set_xticks(range(0,25,2))
f.suptitle('Hormonas, 3 comidas (Fig. 2A): código vs. ecuaciones del texto',fontsize=10); f.tight_layout(); f.savefig('fig1_hormonas.png',dpi=150)
# ---- Fig 2 neurons 24h
names=['AgRP','POMC','PVH','LHA']; qi=[1,2,0,3]
f,ax=plt.subplots(4,2,figsize=(10,7),sharex=True,sharey=True)
x=np.arange(0,24,0.5)+0.25
for c,diet in enumerate(['TRE','CD']):
    for nm in ['publicado','LEPR_POMC_corregido']:
        r=L(diet,24,nm)
        for k,(n,q) in enumerate(zip(names,qi)):
            ax[k,c].plot(x,rate(r,q,24),'-o',ms=2.5,color=COL[nm],lw=1.2,label=LAB[nm] if k==0 else None)
            ax[k,c].set_ylabel(n+'\n(spikes/30 min)',fontsize=8)
    ax[0,c].set_title({'TRE':'3 comidas (Fig. 2A)','CD':'3 comidas + 2 meriendas (Fig. 2B)'}[diet],fontsize=10)
    ax[-1,c].set_xlabel('Tiempo (h)'); ax[-1,c].set_xticks(range(0,25,4))
ax[0,0].legend(fontsize=7); f.tight_layout(); f.savefig('fig2_neuronas_24h.png',dpi=150)
# ---- Fig 3 72h
f,ax=plt.subplots(3,2,figsize=(11,5.5),sharex=True,sharey=True); x=np.arange(0,72,0.5)+0.25
for c,diet in enumerate(['TRE','CD']):
    for nm in ['publicado','LEPR_POMC_corregido']:
        r=L(diet,72,nm)
        for k,(n,q) in enumerate(zip(['POMC','PVH','LHA'],[2,0,3])):
            ax[k,c].plot(x,rate(r,q,72),'-',color=COL[nm],lw=1.1,label=LAB[nm] if k==0 else None); ax[k,c].set_ylabel(n,fontsize=8)
    ax[0,c].set_title({'TRE':'TRE, 72 h (Fig. 4A)','CD':'Dieta común, 72 h (Fig. 4B)'}[diet],fontsize=10)
    ax[-1,c].set_xticks(range(0,73,12)); ax[-1,c].set_xlabel('Tiempo (h)')
ax[0,0].legend(fontsize=7); f.tight_layout(); f.savefig('fig3_72h.png',dpi=150)
# ---- metrics
print('RMSE vs datos (ghrelina, insulina, leptina pM):')
for o in out: print('  %-22s %6.1f %6.1f %6.2f'%o)
def episodes(v,th=10): a=(v>=th).astype(int); return int(((a[1:]-a[:-1])==1).sum()+a[0])
for diet,h in (('TRE',24),('CD',24),('TRE',72),('CD',72)):
    p=L(diet,h,'publicado'); c=L(diet,h,'LEPR_POMC_corregido')
    print(diet,h)
    for n,q in zip(names,qi):
        a,b=rate(p,q,h),rate(c,q,h)
        print('  %-5s total %5d -> %5d (%+.0f%%)  episodios %d -> %d  r=%.2f'%(n,a.sum(),b.sum(),100*(b.sum()-a.sum())/a.sum(),episodes(a),episodes(b),np.corrcoef(a,b)[0,1]))
