import model as M,numpy as np,pickle,sys
V={'publicado':[1,0,0],'LEPR_POMC_corregido':[1,0,1],'leptina_como_texto':[-1,0,0],'ghrelina_como_texto':[1,1,0]}
diet,hours=sys.argv[1],float(sys.argv[2]); names=sys.argv[3].split(',')
for nm in names:
    rec,sp,ns=M.run(M.P,M.meals(diet),np.array(V[nm]),hours*3600.,0.01,1000)
    pickle.dump({'rec':rec,'sp':[sp[q,:ns[q]].copy() for q in range(4)]},open(f'res_{diet}_{int(hours)}_{nm}.pkl','wb'))
    print(diet,hours,nm,ns,flush=True)
