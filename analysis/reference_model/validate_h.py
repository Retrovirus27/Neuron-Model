import pickle,numpy as np,model_h as MH,sys
diet=sys.argv[1]
rec,sp=MH.simulate(diet,24)
ref=pickle.load(open(f'res_{diet}_24_LEPR_POMC_corregido.pkl','rb'))
b=np.arange(0,24.01,0.5)
for i,nm in enumerate(['PVH','AgRP','POMC','LHA']):
    a=np.histogram(ref['sp'][i]/3600,b)[0]; h=np.histogram(sp[i],b)[0]
    print(diet,nm,'total s/h',a.sum(),h.sum(),'max|dif bin|',np.abs(a-h).max())
R=ref['rec']; print(diet,'hormonas max rel dif', np.max(np.abs(R[:,17:22]-rec[:R.shape[0],17:22])/(np.abs(R[:,17:22])+1e-9)))
pickle.dump({'rec':rec,'sp':sp},open(f'resh_{diet}_24.pkl','wb'))
