import sys, numpy as np
HC=219474.6313705; d=np.load(f'ff_{sys.argv[1]}.npz'); w=d['omega']; n=len(w); dQ=float(d['dQ'])
acc={}
for k in range(n):
    D1=(d['Hqp'][k]-d['Hqm'][k])/(2*dQ); D2=(d['Hqp'][k]+d['Hqm'][k]-2*d['Hq0'])/dQ**2
    for i in range(n):
        for j in range(n):
            for key,v in ((tuple(sorted((i,j,k))),D1[i,j]),(tuple(sorted((i,j,k,k))),D2[i,j])):
                acc.setdefault(key,[]).append(v/np.sqrt(np.prod(w[list(key)])))
print('key            mean(cm-1)   spread(cm-1)  [phi in dimensionless q]')
for key,v in sorted(acc.items(), key=lambda t:-abs(np.mean(t[1]))):
    v=np.array(v)*HC
    if abs(v.mean())>2 or np.ptp(v)>1: print(f'{str(key):14s} {v.mean():10.2f} {np.ptp(v):10.2f}  n={len(v)}')
