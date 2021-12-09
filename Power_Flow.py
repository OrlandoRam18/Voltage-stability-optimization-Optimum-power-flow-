# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 19:13:50 2020

@author: orlan_qtpj2q1
"""
def PowerFlowAC(mat1, mat2, mat3,Vminimo,Vmaximo):
    import numpy as np
    from M_YBUS import Y_bus # Función para construir Ybus
    from M_YBUS import mismatch_Power # Función para calcular la diferencia de potencia
    from M_YBUS import Jacobian # Función para calcular la diferencia de potencia
    
    # Funcion de Cal power    
    def cal_power(G,B,nodes,Tetha,V):
        import numpy as np
        P_cal=np.zeros((nodes,1)); Q_cal=np.zeros((nodes,1))
        for k in range(0,nodes):
            for m in range(0,nodes):
                P_cal[k] = P_cal[k] + V[k]*V[m]*( G[k,m]*np.cos(Tetha[k]-Tetha[m]) + B[k,m]*np.sin(Tetha[k]-Tetha[m]))
                Q_cal[k] = Q_cal[k] + V[k]*V[m]*( G[k,m]*np.sin(Tetha[k]-Tetha[m]) - B[k,m]*np.cos(Tetha[k]-Tetha[m]))
        return P_cal, Q_cal
    
    
    # Calcula la matriz de admitancias
    (Y)=Y_bus(mat2,mat1) # Matriz de admitancias
    
    V=np.asarray(mat1[:,1].tolist()); V=V[:, None]# voltajes inciales en pu
    Theta=np.asarray(mat1[:,2].tolist()); Theta=Theta[:, None] # angulo inicial
    #  Generacion y carga de la red
    P_neta = mat1[:,3]-mat1[:,5];  P_neta = P_neta[:, None]
    # Generación y carga reactiva con inyección del shunt
    Q_neta =mat1[:,4]-mat1[:,6]; Q_neta = Q_neta[:, None]
    # Matriz de conductancia y suceptancia
    G=Y.real; B=Y.imag
    # Tipo de bus
    slack=[]; PV=[]; PQ=[]
    for k in range(0,len(mat1)):
        if mat1[k,-1] == 1: slack.append(mat1[k,0])
        if mat1[k,-1] == 2: PV.append(mat1[k,0])
        if mat1[k,-1] == 3: PQ.append(mat1[k,0])
    # Total de buses por tipo
    nslack=len(slack); nPQ=len(PQ); nPV=len(PV)
    tol_max = 10e-6;    # tolerancia de la solucion
    # Metodo iterativo Newton Raphson
    Tol= 1; iter = 1
    while Tol > tol_max:
        [P_cal,Q_cal]=cal_power(G,B,len(mat1),Theta,V); # Esta funcion calcula la potencia
        mismatch_PQ=mismatch_Power(P_neta,Q_neta,P_cal,Q_cal,len(mat1),PQ,nPQ); # Esta funcion calcula la diferencia de potencia
        J=Jacobian(G,B,V,Theta,len(mat1),PQ,nPQ); # Matriz jacobiana de la n-esima iteracion
        mismatch_X=np.linalg.inv(J).dot(mismatch_PQ)
        dTeth=mismatch_X[0:len(mat1)-1]; dV=mismatch_X[len(mat1)-1::]
        # Actualiza V y theta
        Theta[1:len(mat1)]=Theta[1:len(mat1)]+dTeth
        for k in range(0,nPQ):
            n = int(PQ[k]-1)
            V[n]=V[n]+dV[k]
            
        # Check for reactive limits
        if iter >=2 and iter < 7:
            Vma = np.zeros(len(V),dtype=complex)
            for k in range(0,len(V)): Vma[k] = complex(V[k]*np.cos(Theta[k]),V[k]*np.sin(Theta[k]))
            Vma=Vma[:, None]
            Ia= Y.dot(Vma) # Corriente en los buses
            Sa= Vma*Ia.conjugate() # potencia compleja en los nodos            
            for i in range(0,len(mat3)):
                for j in range(0,len(mat1)):
                    if mat1[j,0] == mat3[i,1]:  
                        # Por debajo del limite inferior
                        if Sa[j].imag < mat3[i,5]:
                            V[j] = V[j] + 0.01
                        if V[j] > Vmaximo:
                            V[j] = Vmaximo
                        # Por encima del limite superior
                        if Sa[j].imag > mat3[i,4]:
                            V[j] = V[j] - 0.01
                        if V[j] < Vminimo:
                            V[j] = Vminimo
        iter = iter + 1
        Tol = max(abs(mismatch_PQ))
        
    # Cacula las potencias en los nodos
    Vm=np.zeros(len(V),dtype=complex)
    for k in range(0,len(V)): Vm[k]=complex(V[k]*np.cos(Theta[k]),V[k]*np.sin(Theta[k]))
    Vm=Vm[:, None]
    I= Y.dot(Vm) # Corriente en los buses
    S= Vm*I.conjugate() # potencia compleja en los nodos
    Pgen=[]
    for k in range(0,len(mat1)):
        if mat1[k,-1] ==1: Pgen.append(S[k].real)
        elif mat1[k,-1] ==2: Pgen.append(S[k].real+mat1[k,5])
        elif mat1[k,-1] ==3: Pgen.append(0)
    # Calcula la suceptancia de las lineas
    FromNode = []; ToNode = []
    for i in range(0,len(mat2)):
        FromNode.append(int(mat1[:,0].tolist().index(mat2[i,0])))
        ToNode.append(int(mat1[:,0].tolist().index(mat2[i,1])))
    suceptancia= np.zeros((len(mat2),2),dtype=complex)
    for k in range(0,len(mat2)):
        a=mat2[k,5]
        if a ==0: # en caso de una linea
            b=complex(0,mat2[k,4])
            suceptancia[k,0]=b/2; suceptancia[k,1]=b/2
        else:
            Zpq = complex(mat2[k,2],mat2[k,3]); Ypq=Zpq**-1
            suceptancia[k,0] = (Ypq/a)*((1/a)-1)
            suceptancia[k,1] = Ypq*(1-(1/a))
    # Define admitancia de las lineas
    r= mat2[:,2]; rx= mat2[:,3]; z= np.zeros(len(mat2),dtype=complex)
    for k in range(0,len(mat2)): z[k]=complex(r[k],rx[k])
    y=np.ones(len(mat2))/z
    # Define los flujos de potencia compleja del nodo que envia
    loss_q = []
    Ss=np.zeros(len(mat2),dtype=complex)
    for i in range(0,len(mat2)):
        aux1=Vm[FromNode[i]]; aux2=((Vm[FromNode[i]] - Vm[ToNode[i]])*y[i] + Vm[FromNode[i]]*suceptancia[i,0]).conjugate()
        Ss[i]=aux1*aux2 
        loss_q.append(y[i].imag*(-V[FromNode[i]]**2 - V[ToNode[i]]**2 +2*(Vm[FromNode[i]]*(Vm[ToNode[i]].conjugate())).real))
    # Define los flujos de potencia compleja del nodo receptor
    Sr = np.zeros(len(mat2),dtype = complex)
    for i in range(0,len(mat2)):
        aux3=Vm[ToNode[i]]; aux4=((Vm[ToNode[i]] - Vm[FromNode[i]])*y[i] + Vm[ToNode[i]]*suceptancia[i,1]).conjugate()
        Sr[i]=aux3*aux4 
    # Find the max active power flow         
    Smax = np.zeros(len(mat2))
    for k in range(0,len(mat2)):
        Smax[k] = max(abs(Ss[k]),abs(Sr[k]))

    # Inyecciones activas y reactivas
    Piny = []; Qiny = []; Vgen = []
    for k in range(0,len(mat1)):
        if mat1[k,-1] ==1: 
            Vgen.append(V[k])
            Piny.append(S[k].real+mat1[k,5])
            Qiny.append(S[k].imag+mat1[k,6])
        elif mat1[k,-1] ==2: 
            Vgen.append(V[k])
            Piny.append(S[k].real+mat1[k,5])
            Qiny.append(S[k].imag+mat1[k,6])

    # Define active and reactive power flows
    Pij = np.real(Ss)
    Qij = np.imag(Ss)
    Pji = np.real(Sr)
    Qji = np.imag(Sr)
    return V, Vm, Theta, Pgen, Ss, PQ, PV, slack, Smax, Piny, Qiny, Vgen, Pij, Pji, loss_q


