import pickle,numpy as np,model_h as MH,sys,os
diet=sys.argv[1]; f=f'resh_{diet}_24.pkl'
if os.path.exists(f): d=pickle.load(open(f,'rb')); rec,sp=d['rec'],d['sp']
else: rec,sp=MH.simulate(diet,24); pickle.dump({'rec':rec,'sp':sp},open(f,'wb'))
b=np.arange(0,24.01,0.5); C=[np.histogram(sp[i],b)[0] for i in (1,2,0,3)]
os.makedirs('ref',exist_ok=True)
with open(f'ref/spikes_{diet}_24h.csv','w') as o:
    o.write('t_inicio_h,AgRP,POMC,PVH,LHA\n')
    for k in range(48): o.write('%.1f,%d,%d,%d,%d\n'%(b[k],C[0][k],C[1][k],C[2][k],C[3][k]))
t=np.arange(rec.shape[0])*1000*0.01/3600
with open(f'ref/hormonas_{diet}_24h.csv','w') as o:
    o.write('t_h,GL1,GL2,Ghrelina_pM,Insulina_pM,Leptina_pM\n')
    for h in np.arange(0,24.01,0.5):
        k=int(round(h*3600/10)); k=min(k,rec.shape[0]-1); r=rec[k]
        o.write('%.1f,%.10g,%.10g,%.10g,%.10g,%.10g\n'%(h,r[17],r[18],r[19]/3.3709,r[20],r[21]/2.0043))
print(diet,'ok',[c.sum() for c in C])
