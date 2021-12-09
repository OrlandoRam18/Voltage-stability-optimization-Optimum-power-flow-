# -*- coding: utf-8 -*-
"""
Created on Sun Aug  8 19:13:45 2021

@author: Orlando RB
"""

def Lindex(V,gen,buso,Po,Qo,Line,Sflow,pq,pv,slk,Hshipley,Vcomplejo,Vm,VM):
    import numpy as np
    # Primeramente se calculan las violaciones a los limites operativos
    Pbase = 100; Vbase = 230
    # Revisar violaciones de voltaje
    X = []
    DV=[]; Vmin = Vm; Vmax = VM
    for i in range(0,len(V)):
        if buso[i,-1] == 1 or buso[i,-1] == 2:
            if V[i] < Vmin: DV.append((Vmin-V[i])**2); X.append(1)
            elif V[i] > Vmin and V[i] < Vmax: DV.append(0); X.append(0)
            elif V[i] > Vmax: DV.append((V[i]-Vmax)**2); X.append(1) 
    Dv=sum(DV)*Vbase
    X = sum(X)
    # Revisar violaciones de generación activa
    DG=[]; Y = []
    for i in range(0,len(gen)):
        if Po[i] > gen.iloc[i,2]: DG.append((Pbase*(Po[i]-gen.iloc[i,2]))**2); Y.append(1)
        elif Po[i] > gen.iloc[i,3] and Po[i] < gen.iloc[i,2]: DG.append(0); Y.append(0)
        elif Po[i] < gen.iloc[i,3]: DG.append((Pbase*(gen.iloc[i,3]-Po[i]))**2); Y.append(1)
    Dg=sum(DG)
    Y = sum(Y)

    # Revisar violaciones de generación reactiva
    DQ=[]; Z = []
    for i in range(0,len(gen)):
        if Qo[i] > gen.iloc[i,4]: DQ.append((Pbase*(Qo[i]-gen.iloc[i,4]))**2); Z.append(1)
        elif Qo[i] > gen.iloc[i,5] and Qo[i] < gen.iloc[i,4]: DQ.append(0); Z.append(0)
        elif Qo[i] < gen.iloc[i,5]: DQ.append((Pbase*(gen.iloc[i,5]-Qo[i]))**2); Z.append(1)
    Dq=sum(DQ)
    Z = sum(Z)

    # Revisar violaciones en los limites de transmision
    DPLL=[]; W = []
    for k in range(0,len(Line)):
        if Sflow[k] > Line[k,7]: DPLL.append((Pbase*(Sflow[k]-Line[k,7]))**2); W.append(1)
        else: DPLL.append(0); W.append(0)
    Dl=sum(DPLL)
    W = sum(W)

    
    # Calcula el indice de estabilidad de voltaje de nodos PQ
    Bgen = slk+pv # buses con generación
    lBgen = len(Bgen) # numero de nodos generadores
    lPQ = len(pq)   # Cantidad de nodos PQ
    # Matriz de relación entre buses de carga y nodos PV
    Flg = Hshipley[0:len(pq),len(pq)::]
    
    # vector de voltages de nodos de generación
    Vg = np.zeros(len(Bgen),dtype=complex)
    for j in range(0,len(Bgen)):
        i = buso[:,0].tolist().index(int(Bgen[j]))
        #i = int(Bgen[j]-1)
        Vg[j] = Vcomplejo[i]
    # vector de voltages de nodos de carga
    Vl = np.zeros(len(pq),dtype=complex)
    for j in range(0,len(pq)):
        i = buso[:,0].tolist().index(int(pq[j]))
        #i = int(pq[j]-1)
        Vl[j] = Vcomplejo[i]
    
    # Indice de estabilidad de voltaje
    L = np.zeros(len(pq),dtype=complex)
    for j in range(0,len(pq)):
        for i in range(0,len(Bgen)):
            L[j] = L[j] + (Flg[j,i]*Vg[i])/Vl[j]
        L[j] = abs(1-L[j])
    
    # Define la funcion objetivo    
    Fobj = (Dv+Dg+Dq+Dl) + max(L).real
    return Fobj
