# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 10:18:25 2020

@author: orlan_qtpj2q1
"""


# Construcción de la matriz de admitancias
def Y_bus(line,bus):
    nbuses = len(bus); 
    import numpy as np
    Y=np.zeros((nbuses,nbuses),dtype=complex)
    orden=[]
    for k in range(1,nbuses+1): orden.append(k)
    for p in range(0,len(line)):
        BusP = bus[:,0].tolist().index(line[p,0]) 
        BusQ = bus[:,0].tolist().index(line[p,1])
        a=line[p,5]; # Tap value for the  p iteration
        if a > 0: 
            yl=(1/(complex(line[p,2],line[p,3]))); # admitancia de la linea
            Ad=complex(0,line[p,4]/2) # line charging
            Y[BusP,BusQ]=Y[BusP,BusQ]-yl/a; # elemento no diagonal
            Y[BusQ,BusP]=Y[BusP,BusQ]; # se declara la simetria en elementos fuera de la diagonal
            Y[BusP,BusP]=Y[BusP,BusP]+(yl/a)+((1/a)*(1/a-1)*yl)+Ad; # Equivalente de admitancia en el bus P + line charging
            Y[BusQ,BusQ]=Y[BusQ,BusQ]+(yl/a)+(1-1/a)*yl+Ad; # Equivalente de admitancia en el bus Q + line charging
        else: # En caso de una línea
            yl=(1/(complex(line[p,2],line[p,3]))); # admitancia de la linea
            Ad=complex(0,line[p,4]/2) # line charging
            Y[BusP,BusQ]=Y[BusP,BusQ]-yl # elemento no diagonal
            Y[BusQ,BusP]=Y[BusP,BusQ]; # se declara la simetria en elementos fuera de la diagonal
            Y[BusP,BusP]=Y[BusP,BusP]+yl # elemento diagonal
            Y[BusQ,BusQ]=Y[BusQ,BusQ]+yl # elemento diagonal
            c=line[p,4]; # line charging de toda la linea
            if c > 0:
                Y[BusP,BusP]= Y[BusP,BusP]+Ad; # Añade el valor del charging al elemento diagonal
                Y[BusQ,BusQ]= Y[BusQ,BusQ]+Ad; # Añade el valor del charging al elemento diagonal
    for q in range(0,nbuses):
        Y[int(bus[q,0]-1),int(bus[q,0]-1)] += complex(bus[q,7],bus[q,8])
    return Y


# Construccion de la matriz de admitancias modificada
def YBUS_mod(line,bus,pq,pv,slk):
    nbuses = len(bus)
    import numpy as np
    Y=np.zeros((nbuses,nbuses),dtype=complex)
    Bgen = slk+pv
    indice_or = pq+Bgen
    for p in range(0,len(line)):
        BusP = indice_or.index(line[p,0]) 
        BusQ = indice_or.index(line[p,1])
        a=line[p,5]; # Tap value for the  p iteration
        if a > 0: 
            yl=(1/(complex(line[p,2],line[p,3]))); # admitancia de la linea
            Ad=complex(0,line[p,4]/2) # line charging
            Y[BusP,BusQ]=Y[BusP,BusQ]-yl/a; # elemento no diagonal
            Y[BusQ,BusP]=Y[BusP,BusQ]; # se declara la simetria en elementos fuera de la diagonal
            Y[BusP,BusP]=Y[BusP,BusP]+(yl/a)+((1/a)*(1/a-1)*yl)+Ad; # Equivalente de admitancia en el bus P + line charging
            Y[BusQ,BusQ]=Y[BusQ,BusQ]+(yl/a)+(1-1/a)*yl+Ad; # Equivalente de admitancia en el bus Q + line charging
        else: # En caso de una línea
            yl=(1/(complex(line[p,2],line[p,3]))); # admitancia de la linea
            Ad=complex(0,line[p,4]/2) # line charging
            Y[BusP,BusQ]=Y[BusP,BusQ]-yl # elemento no diagonal
            Y[BusQ,BusP]=Y[BusP,BusQ]; # se declara la simetria en elementos fuera de la diagonal
            Y[BusP,BusP]=Y[BusP,BusP]+yl # elemento diagonal
            Y[BusQ,BusQ]=Y[BusQ,BusQ]+yl # elemento diagonal
            c=line[p,4]; # line charging de toda la linea
            if c > 0:
                Y[BusP,BusP]= Y[BusP,BusP]+Ad; # Añade el valor del charging al elemento diagonal
                Y[BusQ,BusQ]= Y[BusQ,BusQ]+Ad; # Añade el valor del charging al elemento diagonal
    for q in range(0,nbuses):
        BusP = indice_or.index(bus[q,0])
        Y[BusP,BusP] += complex(bus[BusP,7],bus[BusP,8])
    return Y,indice_or

def Y_busDC(line,bus,nbuses):
    import numpy as np
    Y=np.zeros((nbuses,nbuses))
    orden=[]
    line[:,2] = np.zeros(len(line)) # Ignora la parte real de la linea
    for k in range(1,nbuses+1): orden.append(k)
    for p in range(0,len(line)):
        BusP = bus[:,0].tolist().index(line[p,0]) 
        BusQ = bus[:,0].tolist().index(line[p,1])
        a=line[p,5]; # Tap value for the  p iteration
        if a > 0: 
            yl=(1/((line[p,2]+line[p,3]))); # admitancia de la linea
            Ad=(line[p,4]/2) # line charging
            Y[BusP,BusQ]=Y[BusP,BusQ]-yl/a; # elemento no diagonal
            Y[BusQ,BusP]=Y[BusP,BusQ]; # se declara la simetria en elementos fuera de la diagonal
            Y[BusP,BusP]=Y[BusP,BusP]+(yl/a)+((1/a)*(1/a-1)*yl)+Ad; # Equivalente de admitancia en el bus P + line charging
            Y[BusQ,BusQ]=Y[BusQ,BusQ]+(yl/a)+(1-1/a)*yl+Ad; # Equivalente de admitancia en el bus Q + line charging
        else: # En caso de una línea
            yl=(1/((line[p,2]+line[p,3]))); # admitancia de la linea
            Ad=(line[p,4]/2) # line charging
            Y[BusP,BusQ]=Y[BusP,BusQ]-yl # elemento no diagonal
            Y[BusQ,BusP]=Y[BusP,BusQ]; # se declara la simetria en elementos fuera de la diagonal
            Y[BusP,BusP]=Y[BusP,BusP]+yl # elemento diagonal
            Y[BusQ,BusQ]=Y[BusQ,BusQ]+yl # elemento diagonal
            c=line[p,4]; # line charging de toda la linea
            if c > 0:
                Y[BusP,BusP]= Y[BusP,BusP]+Ad; # Añade el valor del charging al elemento diagonal
                Y[BusQ,BusQ]= Y[BusQ,BusQ]+Ad; # Añade el valor del charging al elemento diagonal
    return Y
    



def mismatch_Power(P_neta,Q_neta,P_cal,Q_cal,nodes,PQ,nPQ):
    import numpy as np
    dP = P_neta-P_cal; dQ = Q_neta-Q_cal; l=0
    mismatch_Q=[]
    for k in range(0,nPQ):
        n = int(PQ[k]-1)
        #n = PQ.index(PQ[k])
        #n=int(PQ[k]-1); 
        mismatch_Q.append(dQ[n]); l=l+1
    mismatch_P = 1*dP[1:nodes]; # vector de error de potencia sin el bus slack
    mismatch_PQ = np.concatenate((mismatch_P, np.asarray(mismatch_Q))) # Mismatch vector  [dP;dQ]
    return mismatch_PQ




def Jacobian(G,B,V,Tetha,nodes,PQ,nPQ):
    import numpy as np
    J1 = np.zeros((nodes-1,nodes-1));
    J2 = np.zeros((nodes-1,nPQ));
    J3 = np.zeros((nPQ,nodes-1));
    J4 = np.zeros((nPQ,nPQ));
    # J1 [dPk/dTethaj]
    for k in range(1,nodes):
        for j in range(1,nodes):
            if j == k: # elemento de la diagonal
                for m in range(0,nodes): # suma de cada elemento
                    J1[k-1,j-1]=J1[k-1,j-1]+(V[k]*V[m])*(-G[k,m]*np.sin(Tetha[k]-Tetha[m])+B[k,m]*np.cos(Tetha[k]-Tetha[m]))
                J1[k-1,j-1]=J1[k-1,j-1]-(V[k]**2)*B[k,k] # sustrae el elemento cuando m= k, dado que su derivada es 0
            else: # para elementos no diagonales
                J1[k-1,j-1]=V[k]*V[j]*(G[k,j]*np.sin(Tetha[k]-Tetha[j])-B[k,j]*np.cos(Tetha[k]-Tetha[j]))
    ##
    # J2 [dPk/dV_j]
    for k in range(1,nodes):
        for j in range(0,nPQ):
            m = int(PQ[j]-1)
            #m = PQ.index(PQ[j])
            if m == k: # elemento de la diagonal
                for m in range(0,nodes):
                    J2[k-1,j]=J2[k-1,j]+V[m]*(G[k,m]*np.cos(Tetha[k]-Tetha[m])+B[k,m]*np.sin(Tetha[k]-Tetha[m]))
                J2[k-1,j]=J2[k-1,j]+(V[k])*G[k,k]
            else: # elemento no diagonal
                J2[k-1,j]=V[k]*(G[k,m]*np.cos(Tetha[k]-Tetha[m])+B[k,m]*np.sin(Tetha[k]-Tetha[m]))
    ##
    # J3 [dQk/dTeth_j]
    for k in range(0,nPQ):
        #n = PQ.index(PQ[k])
        n = int(PQ[k]-1)
        for j in range(1,nodes):
            if j == n:
                for m in range(0,nodes):
                    J3[k,j-1]=J3[k,j-1]+(V[n]*V[m])*(G[n,m]*np.cos(Tetha[n]-Tetha[m])+B[n,m]*np.sin(Tetha[n]-Tetha[m]))
                J3[k,j-1]=J3[k,j-1]-(V[n]**2)*G[n,n]
            else: 
                J3[k,j-1]=V[n]*V[j]*(-G[n,j]*np.cos(Tetha[n]-Tetha[j])-B[n,j]*np.sin(Tetha[n]-Tetha[j]))
    ##
    # J4 [dQk/dV_j]
    for k in range(0,nPQ):
        #n = PQ.index(PQ[k])
        n = int(PQ[k]-1)
        for j in range(0,nPQ):
            m = int(PQ[j]-1)
            #m = PQ.index(PQ[j])
            if m == n:
                for m in range(0,nodes):
                    J4[k,j] = J4[k,j] + V[m]*(G[n,m]*np.sin(Tetha[n]-Tetha[m]) - B[n,m]*np.cos(Tetha[n]-Tetha[m]))
                J4[k,j]=J4[k,j]-V[n]*B[n,n]
            else:
                J4[k,j]=V[n]*(G[n,m]*np.sin(Tetha[n]-Tetha[m])- B[n,m]*np.cos(Tetha[n]-Tetha[m]))
    J=np.concatenate((np.hstack((J1,J2)),np.hstack((J3,J4))))
    return J
