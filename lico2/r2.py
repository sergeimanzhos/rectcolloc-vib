import numpy as np
def r2(idx,d):
    g=2.0
    for _ in range(100): g=(1+g)**(1/(d+1))
    a=(1/g)**np.arange(1,d+1); return (0.5+np.outer(np.asarray(idx,float),a))%1.0
