# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 19:30:34 2021

@author: Orlando RB
"""

# Unidad de recursos de generación
# Codigo principal --- Optimización de estabilidad de perfiles de voltaje


# Importa librerias necesarias
import pandas as pd; import numpy as np; import random; from datetime import date
import time

# Leectura de la base de datos
print("")
today = date.today().strftime("%d/%m/%Y")
print(date.today() ,'  Inicia lectura del sistema a analizar')
print('             ',date.today() ,'  ',time.strftime("%X"))
Data_collection = pd.read_excel("14 bus data.xlsx", sheet_name=None,  header=0)
print(date.today() ,'  Finaliza lectura del sistema a analizar')

# Define las 4 matrices de datos
LINE = 1*Data_collection['Line']; BUS = 1*Data_collection['Bus']; 
GEN = 1*Data_collection['Gen']; SHUNT = 1*Data_collection['Shunts']; 


# # ++++++++++++++++ Evalua la red original y determina el plan de expansión ++++++++++++++++
print("")
# Inicia la población del algoritmo de optimización
poblacion = 10
if poblacion*(len(GEN)+len(SHUNT)) <= 60: pop = 60
else: pop = 120
# Define la población inicial de los voltajes de los nodos PV
print('   Población inicial de los voltajes de los nodos PV')
Vmin = 0.93; Vmax = 1.07
x_vgen = np.zeros((pop,len(GEN))) 
for k in range(0,pop):
    for j in range(0,len(GEN)):
        x_vgen[k,j] = Vmin + (Vmax - Vmin)*np.random.rand(1)
print('   Propuesta de los niveles de generación')
x_gen = np.zeros((pop,len(GEN)))
for k in range(0,pop):
    for j in range(0,len(GEN)):
        x_gen[k,j] = GEN.iloc[j,3] + (GEN.iloc[j,2] - GEN.iloc[j,3])*np.random.rand(1)
print('   Corrigiendo generación')
for k in range(0,pop):
    # si la generación es mayor que la carga
    if sum(x_gen[k,:]) > 1.01*sum(BUS.iloc[:,5]):
        alfa = sum(BUS.iloc[:,5])/sum(x_gen[k,:])
        x_gen[k,:] = alfa*x_gen[k,:]
    # si la carga es mayor que la generación
    elif sum(BUS.iloc[:,5]) > sum(x_gen[k,:]):
        beta = sum(x_gen[k,:])/sum(BUS.iloc[:,5])
        x_gen[k,:] = x_gen[k,:]/beta
                
# Define la población inicial de los Taps de los transformadores
# Busca transformadores en el sistema
TAP = []; a = LINE.columns.tolist().index('Tipo de línea')
for k in range(0,len(LINE)):
    if LINE.iloc[k,a] > 0: TAP.append(LINE.iloc[k,::])
TAP = pd.DataFrame(TAP)
if len(TAP) > 0:
    print('   Población inicial de las posiciones de los taps')
    Tmin = 0.9; Tmax = 1.1
    x_tap = np.zeros((pop,len(TAP)))
    for k in range(0,pop):
        for j in range(0,len(TAP)):
            x_tap[k,j] = (random.randint(int(Tmin*100),int(Tmax*100))/100)
else:
    print('   No se encontraron transformadores en la red')
# # Define la población inicial de los shunts
if len(SHUNT) > 0:
    print('   Población inicial de las suceptancias shunt')
    Smax = np.zeros(len(SHUNT))
    Smin = np.zeros(len(SHUNT))
    x_shunt = np.zeros((pop,len(SHUNT)))
    maxind = SHUNT.columns.tolist().index('Bsvc max pu')
    minind = SHUNT.columns.tolist().index('Bsvc min pu')
    for k in range(0,len(SHUNT)):
        Smax[k] = SHUNT.iloc[k,maxind]
        Smin[k] = SHUNT.iloc[k,minind]
    for k in range(0,pop):
        for j in range(0,len(SHUNT)):
            x_shunt[k,j] = Smin[j] + (Smax[j] - Smin[j])*np.random.rand(1)
else:
    print('   No se encontraron elementos shunt en la red')
# Total de iteraciones
n_itera = 10
# Función para determinar elnumero de subconjuntos
def Mvalue(m,Obf1,Obf0):
    if min(Obf0) <= min(Obf1): m = m+1
    elif min(Obf0) > min(Obf1) and m >=3: m = m - 1
    # Corrige valor de m en caso de ser muy grande
    if m > 6: m = 6
    return m
# Numero inicial de subconjuntos
M = 2

# Matrices de almacenamiento de información
bsof=[]
xigen = np.zeros((n_itera,len(GEN))); xivgen = np.zeros((n_itera,len(GEN)))
Obfun=np.zeros((n_itera,pop)); Obfunnew=np.zeros((n_itera,pop))
XPG=np.zeros((n_itera,pop,len(GEN))) # guarda datos de cadidatos de generación
XVG=np.zeros((n_itera,pop,len(GEN))) # guarda datos de voltajes candidatos
if len (TAP) > 0:
    XTAP = np.zeros((n_itera,pop,len(TAP))) # guarda datos de TAPS candidatos
    xitap = np.zeros((n_itera,len(TAP)))
if len (SHUNT) > 0:
    XSHUNT = np.zeros((n_itera,pop,len(SHUNT))) # guarda datos de suceptancias candidatos
    xishunt = np.zeros((n_itera,len(SHUNT)))

# Llama a las funciones del proceso iterativo
from Power_Flow import PowerFlowAC # Función para calculo de flujos AC
from M_YBUS import YBUS_mod # Función para construir Ybus modificada
from Invparcial import shipley # Función para construir inversa parcial
from Fobjetivo import Lindex # Función para construir inversa parcial
from ActualizaX import UpdateDesitionVariablesAA # Función para actualizar 2 variables de desición
from ActualizaX import UpdateDesitionVariablesBB # Función para actualizar 3 variables de desición
from ActualizaX import UpdateDesitionVariablesCC # Función para actualizar 4 variables de desición


print("")
print(date.today() ,'  Inicia proceso iterativo')
print('             ',date.today() ,'  ',time.strftime("%X"))
# +++++++++++++++++++++++++++ Inicia proceso iterativo +++++++++++++++++++++++
for k in range(0,n_itera):
    
    # Crea los conjuntos de pobalciones
    m=[]
    for n in range(0,M):
        if n==0: m.append(round(pop/M))
        #elif n == int(M-1): m.append(int(pop))
        else:  m.append(round((n+1)*pop/M))
    
    for p in range(0,pop):
        # Matrices auxiliares de datos
        line = 1*LINE; bus = 1*BUS
        shunt = 1*SHUNT; Gen = 1*GEN
        # asigna columna de generacion a cero
        bus.iloc[:,3] = 1*np.zeros(len(BUS))
        # asigna columna de dispositivos de compensación a cero
        bus.iloc[:,8] = 1*np.zeros(len(BUS))
        
        # Modifica la matriz line de acuerdo a la poblacion de los TAPS
        if len(TAP) > 0:
            for i in range(0,len(TAP)):
                for j in range(0,len(line)):
                    if TAP.iloc[i,0] == line.iloc[j,0] and TAP.iloc[i,1] == line.iloc[j,1]:
                        # Modifica el valor del tap
                        line.iloc[j,a] = x_tap[p,i]
        # Modifica la matriz bus de acuerdo a los dispositivos shunt
        if len(SHUNT) > 0:
            for i in range(0,len(SHUNT)):
                for j in range(0,len(bus)):
                    if bus.iloc[j,0] == shunt.iloc[i,1]:
                        # Modifica el valor de la suceptancia
                        # solo en esta ocasión se utiliza la solución inicial del voltage
                        bus.iloc[j,8] = x_shunt[p,i]
        # Modifica la matriz bus de acuerdo a la población de voltajes de generación
        for i in range(0,len(GEN)):
            for j in range(0,len(bus)):
                if bus.iloc[j,0] == GEN.iloc[i,1]:
                    # Modifica el valor del voltaje
                    bus.iloc[j,1] = x_vgen[p,i]
        # Modifica la matriz bus de acuerdo a la población de generación
        for i in range(0,len(GEN)):
            for j in range(0,len(bus)):
                if bus.iloc[j,0] == GEN.iloc[i,1]:
                    # Modifica el valor del generación
                    bus.iloc[j,3] += x_gen[p,i]
                    
        busa=1*np.asarray(bus); linea=1*np.asarray(line); Gena = 1*np.asarray(Gen) # matrices auxiliares
        # Calculo de la solucion de flujos AC
        [Vm,Vcomp,Theta_n,_,Ssend,PQ,PV,Slack,Sflujo,PI,QI,Vgn,_,_,_] = PowerFlowAC(busa,linea,Gena,Vmin,Vmax)
        # Reorganizando la matriz de admitancias
        [Ymod,indexor] = YBUS_mod(linea,busa,PQ,PV,Slack)
        # Inversión parcial de la matriz de admitancias modificada
        H = 1*Ymod
        for Piv in range(0,len(PQ)):
            H = shipley(H,Piv)
        # Calculo de la función objettivo
        Fa = Lindex(Vm,Gen,busa,PI,QI,linea,Sflujo,PQ,PV,Slack,H,Vcomp,Vmin,Vmax)
        # almacena valores de la función objetivo
        Obfun[k,p] = Fa
        
    # Define los nuevos valores de las variables de desicion
    for n in range(0,M):
        if n == 0:
            if len(TAP) > 0 and len (SHUNT) > 0:
                [x2a,x3a,x6a,x7a] = UpdateDesitionVariablesCC(Obfun[k,0:m[0]],x_gen[0:m[0],:],x_vgen[0:m[0],:],x_tap[0:m[0],:],x_shunt[0:m[0],:])
            elif len(TAP) > 0 and len (SHUNT) == 0:
                [x2a,x3a,x6a] = UpdateDesitionVariablesBB(Obfun[k,0:m[0]],x_gen[0:m[0],:],x_vgen[0:m[0],:],x_tap[0:m[0],:])
            elif len(TAP) == 0 and len (SHUNT) > 0:
                [x2a,x3a,x7a] = UpdateDesitionVariablesBB(Obfun[k,0:m[0]],x_gen[0:m[0],:],x_vgen[0:m[0],:],x_shunt[0:m[0],:])
            elif len(TAP) == 0 and len (SHUNT) == 0:
                [x2a,x3a] = UpdateDesitionVariablesAA(Obfun[k,0:m[0]],x_gen[0:m[0],:],x_vgen[0:m[0],:])
        else:
            if len(TAP) > 0 and len (SHUNT) > 0:
                [x2a,X3a,x6a,x7a] = UpdateDesitionVariablesCC(Obfun[k,m[n-1]:m[n]],x_gen[m[n-1]:m[n],:],x_vgen[m[n-1]:m[n],:],x_tap[m[n-1]:m[n],:],x_shunt[m[n-1]:m[n],:])
            elif len(TAP) > 0 and len (SHUNT) == 0:
                [x2a,X3a,x6a] = UpdateDesitionVariablesBB(Obfun[k,m[n-1]:m[n]],x_gen[m[n-1]:m[n],:],x_vgen[m[n-1]:m[n],:],x_tap[m[n-1]:m[n],:])
            elif len(TAP) == 0 and len (SHUNT) > 0:
                [x2a,X3a,x7a] = UpdateDesitionVariablesBB(Obfun[k,m[n-1]:m[n]],x_gen[m[n-1]:m[n],:],x_vgen[m[n-1]:m[n],:],x_shunt[m[n-1]:m[n],:])
            elif len(TAP) == 0 and len (SHUNT) == 0:
                [x2a,X3a] = UpdateDesitionVariablesAA(Obfun[k,m[n-1]:m[n]],x_gen[m[n-1]:m[n],:],x_vgen[m[n-1]:m[n],:])
        if n ==0: 
            if len(TAP) > 0 and len (SHUNT) > 0:
                x2, x3 ,x6, x7 = x2a, x3a, x6a, x7a
            elif len(TAP) > 0 and len (SHUNT) == 0:
                x2, x3 ,x6 = x2a, x3a, x6a
            elif len(TAP) == 0 and len (SHUNT) > 0:
                x2, x3 ,x7 = x2a, x3a, x7a
            elif len(TAP) == 0 and len (SHUNT) == 0:
                x2, x3 = x2a, x3a
        else:
            if len(TAP) > 0 and len (SHUNT) > 0:
                x2 = np.vstack([x2a, x2]); x3 = np.vstack([x3a, x3]); x6 = np.vstack([x6a, x6]); x7 = np.vstack([x7a, x7])
            elif len(TAP) > 0 and len (SHUNT) == 0:
                x2 = np.vstack([x2a, x2]); x3 = np.vstack([x3a, x3]); x6 = np.vstack([x6a, x6])
            elif len(TAP) == 0 and len (SHUNT) > 0:
                x2 = np.vstack([x2a, x2]); x3 = np.vstack([x3a, x3]); x7 = np.vstack([x7a, x7])
            elif len(TAP) == 0 and len (SHUNT) == 0:
                x2 = np.vstack([x2a, x2]); x3 = np.vstack([x3a, x3])
    
    # Corrige todas las particulas que esten sobrepasando sus valores limite
    # Potencia de los generadores
    for jj in range(0,len(Gen)):
        for i in range(0,pop):
            if x2[i,jj] < GEN.iloc[jj,3]: x2[i,jj] = 1*Gena[jj,3]
            if x2[i,jj] > GEN.iloc[jj,2]: x2[i,jj] = 1*Gena[jj,3]        
    x_gen = 1*x2    
    # corrige la potencia de los generadores
    for ka in range(0,pop):
        # si la generación es mayor que la carga
        if sum(x_gen[ka,:]) > 1.01*sum(BUS.iloc[:,5]):
            alfa = sum(BUS.iloc[:,5])/sum(x_gen[ka,:])
            x_gen[ka,:] = alfa*x_gen[ka,:]
        # # si la carga es mayor que la generación
        # elif sum(BUS.iloc[:,5]) > sum(x_gen[ka,:]):
        #     beta = sum(x_gen[ka,:])/sum(BUS.iloc[:,5])
        #     x_gen[ka,:] = x_gen[ka,:]*(1/beta)
    
    # Voltajes de los generadores
    for jj in range(0,len(Gen)):
        for i in range(0,pop):
            if x3[i,jj] < Vmin: x3[i,jj] = 1*Vmin
            if x3[i,jj] > Vmax: x3[i,jj] = 1*Vmax
    x_vgen = 1*x3
    # Valores de los taps
    if len(TAP) > 0:
        # Se convierte a variable discreta
        for jj in range(0,len(TAP)):
            for i in range(0,pop):
                x6[i,jj] = round(x6[i,jj]*100)
                x6[i,jj] = x6[i,jj]/100
        for jj in range(0,len(TAP)):
            for i in range(0,pop):
                if x6[i,jj] < Tmin: x6[i,jj] = 1*Tmin
                if x6[i,jj] > Tmax: x6[i,jj] = 1*Tmax            
        x_tap = 1*x6
    # Valores de las suceptacnias de los shunt
    if len(SHUNT) > 0:
        for jj in range(0,len(SHUNT)):
            for i in range(0,pop):
                if x7[i,jj] < Smin[jj]: x7[i,jj] = 1*Smin[jj]
                if x7[i,jj] > Smax[jj]: x7[i,jj] = 1*Smax[jj]
        x_shunt = 1*x7

    for p in range(0,pop):     
        # Matrices auxiliares de datos
        line = 1*LINE; bus = 1*BUS
        shunt = 1*SHUNT; Gen = 1*GEN
        # asigna columna de generacion a cero
        bus.iloc[:,3] = 1*np.zeros(len(BUS))
        # asigna columna de dispositivos de compensación a cero
        bus.iloc[:,8] = 1*np.zeros(len(BUS))
        
        # Modifica la matriz line de acuerdo a la poblacion de los TAPS
        if len(TAP) > 0:
            for i in range(0,len(TAP)):
                for j in range(0,len(line)):
                    if TAP.iloc[i,0] == line.iloc[j,0] and TAP.iloc[i,1] == line.iloc[j,1]:
                        # Modifica el valor del tap
                        line.iloc[j,a] = 1*x_tap[p,i]
        # Modifica la matriz bus de acuerdo a los dispositivos shunt
        if len(SHUNT) > 0:
            for i in range(0,len(SHUNT)):
                for j in range(0,len(bus)):
                    if bus.iloc[j,0] == shunt.iloc[i,1]:
                        # Modifica el valor de la suceptancia
                        # solo en esta ocasión se utiliza la solución inicial del voltage
                        bus.iloc[j,8] = 1*x_shunt[p,i]
        # Modifica la matriz bus de acuerdo a la población de voltajes de generación
        for i in range(0,len(GEN)):
            for j in range(0,len(bus)):
                if bus.iloc[j,0] == GEN.iloc[i,1]:
                    # Modifica el valor del voltaje
                    bus.iloc[j,1] = 1*x_vgen[p,i]
        # Modifica la matriz bus de acuerdo a la población de generación
        for i in range(0,len(GEN)):
            for j in range(0,len(bus)):
                if bus.iloc[j,0] == GEN.iloc[i,1]:
                    # Modifica el valor del generación
                    bus.iloc[j,3] += 1*x_gen[p,i]
                    
        busb=1*np.asarray(bus); lineb=1*np.asarray(line); Genb = 1*np.asarray(Gen) # matrices auxiliares
        # Calculo de la solucion de flujos AC
        [Vmag,Vcomp,Theta_n,_,Ssend,PQ,PV,Slack,Sflujo,Pin,Qin,Vgn,_,_,_] = PowerFlowAC(busb,lineb,Genb,Vmin,Vmax)
        # Almacena voltajes de resultado de flujos optimos
        #x_vgen[p,:] = Vgn
        # Reorganizando la matriz de admitancias
        [Ymodi,indexor] = YBUS_mod(lineb,busb,PQ,PV,Slack)
        # Inversión parcial de la matriz de admitancias modificada
        Ha = 1*Ymodi
        for Pivv in range(0,len(PQ)):
            Ha = shipley(Ha,Pivv)
        # Calculo de la función objettivo
        Fb = Lindex(Vmag,Gen,busb,Pin,Qin,lineb,Sflujo,PQ,PV,Slack,Ha,Vcomp,Vmin,Vmax)
        # almacena valores de la función objetivo
        Obfunnew[k,p] = Fb
    # Guarda los valores de las variables de desicion para cada iteracion
    XPG[k,:,:]=1*x_gen;   XVG[k,:,:]=1*x_vgen        
    if len(TAP) > 0:
        XTAP[k,:,:] = 1*x_tap
    if len(SHUNT) > 0:
        XSHUNT[k,:,:] = 1*x_shunt
    if k > 0:
        for i in range(0,pop):
            if Obfunnew[k-1,i] < Obfunnew[k,i]:
                x_gen[i,:] = 1*XPG[k-1,i,:]; x_vgen[i,:] = 1*XVG[k-1,i,:]
                if len(TAP) > 0:
                    x_tap[i,:] = 1*XTAP[k-1,i,:]
                if len(SHUNT) > 0:
                    x_shunt[i,:] = 1*XSHUNT[k-1,i,:]
                Obfunnew[k,i] = Obfunnew[k-1,i]
                # En caso de que necesitemos los valores de iteraciones previas
                # Deberemos cambiat nuestras matrices de almacenamiento
                XPG[k,i,:] = 1*x_gen[i,:]; XVG[k,i,:] = 1*x_vgen[i,:]
                if len(TAP) > 0:
                    XTAP[k,i,:] = 1*x_tap[i,:]
                if len(SHUNT) > 0:
                    XSHUNT[k,i,:] = 1*x_shunt[i,:]
    # best solution at each iteration
    bsof.append(min(Obfunnew[k,:]))    
    print(k+1,"   ","Lindex: ", "{:.8f}".format(bsof[k]),'  ',date.today() ,'  ',time.strftime("%X"))
    # Encuentra los valores de las variables de desicion asociados a la mejor solucion
    for i in range(0,pop):
        if bsof[k] == Obfunnew[k,i]:
            xigen[k,:] = 1*x_gen[i,:]; xivgen[k,:] = 1*x_vgen[i,:]
            if len(TAP) > 0:
                xitap[k,:] = 1*x_tap[i,:]
            if len(SHUNT) > 0:
                xishunt[k,:] = 1*x_shunt[i,:]


    # Realiza el ajuste con la poblacion elite
    soluciones=[]
    for n in range(0,M):
        if n ==0: soluciones.append(min(Obfunnew[k,0:m[0]]))
        else:  soluciones.append(min(Obfunnew[k,m[n-1]:m[n]]))
    # Indices de las soluciones de los subconjuntos M
    indice_min = np.asarray(soluciones).argmin()
    indice_max = np.asarray(soluciones).argmax()
    # AJUSTE CON LA POBLACION ELITE
    if indice_min == 0 and indice_max > 0:
        # Modifica las variables locales
        x_gen[m[indice_max-1]:m[indice_max],:] = x_gen[0:m[0],:]
        x_vgen[m[indice_max-1]:m[indice_max],:] = x_vgen[0:m[0],:]
        if len(TAP) > 0:
            x_tap[m[indice_max-1]:m[indice_max],:] = x_tap[0:m[0],:]
        if len(SHUNT) > 0:
            x_shunt[m[indice_max-1]:m[indice_max],:] = x_shunt[0:m[0],:]
    elif indice_min > 0 and indice_max == 0:
        x_gen[0:m[0],:] = x_gen[m[indice_min-1]:m[indice_min],:]
        x_vgen[0:m[0],:] = x_vgen[m[indice_min-1]:m[indice_min],:]
        if len(TAP) > 0:
            x_tap[0:m[0],:] = x_tap[m[indice_min-1]:m[indice_min],:]
        if len(SHUNT) > 0:
            x_shunt[0:m[0],:] = x_shunt[m[indice_min-1]:m[indice_min],:]
    elif indice_min > 0 and indice_max > 0:
        x_gen[m[indice_max-1]:m[indice_max],:] = x_gen[m[indice_min-1]:m[indice_min],:]
        x_vgen[m[indice_max-1]:m[indice_max],:] = x_vgen[m[indice_min-1]:m[indice_min],:]
        if len(TAP) > 0:
            x_tap[m[indice_max-1]:m[indice_max],:] = x_tap[m[indice_min-1]:m[indice_min],:]
        if len(SHUNT) > 0:
            x_shunt[m[indice_max-1]:m[indice_max],:] = x_shunt[m[indice_min-1]:m[indice_min],:]

    # Define un nuevo valor de m
    if k > 0:
        M = Mvalue(M,Obfunnew[k,:],Obfunnew[k-1,:])

print("")
print(date.today() ,'  Finaliza proceso iterativo')
print('             ',date.today() ,'  ',time.strftime("%X"))  
# Matrices auxiliares de datos
line = 1*LINE; bus = 1*BUS
shunt = 1*SHUNT; Gen = 1*GEN
# asigna columna de generacion a cero
bus.iloc[:,3] = 1*np.zeros(len(BUS))
# asigna columna de dispositivos de compensación a cero
bus.iloc[:,8] = 1*np.zeros(len(BUS))

# Con nuevos valores optimos modifica bus y line
# Modifica la matriz line de acuerdo a la poblacion de los TAPS
if len(TAP) > 0:
    for i in range(0,len(TAP)):
        for j in range(0,len(line)):
            if TAP.iloc[i,0] == line.iloc[j,0] and TAP.iloc[i,1] == line.iloc[j,1]:
                # Modifica el valor del tap
                line.iloc[j,a] = 1*xitap[-1,i]
# Modifica la matriz bus de acuerdo a los dispositivos shunt
if len(SHUNT) > 0:
    for i in range(0,len(SHUNT)):
        for j in range(0,len(bus)):
            if bus.iloc[j,0] == shunt.iloc[i,1]:
                # Modifica el valor de la suceptancia
                # solo en esta ocasión se utiliza la solución inicial del voltage
                bus.iloc[j,8] = 1*xishunt[-1,i]
# Modifica la matriz bus de acuerdo a la población de voltajes de generación
for i in range(0,len(GEN)):
    for j in range(0,len(bus)):
        if bus.iloc[j,0] == GEN.iloc[i,1]:
            # Modifica el valor del voltaje
            bus.iloc[j,1] = 1*xivgen[-1,i]
# Modifica la matriz bus de acuerdo a la población de generación
for i in range(0,len(GEN)):
    for j in range(0,len(bus)):
        if bus.iloc[j,0] == GEN.iloc[i,1]:
            # Modifica el valor del generación
            bus.iloc[j,3] += 1*xigen[-1,i]
            
busc=1*np.asarray(bus); linec=1*np.asarray(line); Genc = 1*np.asarray(Gen) # matrices auxiliares
# Calculo de la solucion de flujos AC
[Vmag,Vcomp,Theta_n,_,Ssend,PQ,PV,Slack,Sflujo,Pi,Qi,Vgn,_,_,_] = PowerFlowAC(busc,linec,Genc,Vmin,Vmax)
# Reorganizando la matriz de admitancias
[Ymodx,indexor] = YBUS_mod(linec,busc,PQ,PV,Slack)
# Inversión parcial de la matriz de admitancias modificada
Haa = 1*Ymodx
for Pivv in range(0,len(PQ)):
    Haa = shipley(Haa,Pivv)
# función final
# Calculo de la función objettivo
Fc = Lindex(Vmag,Gen,busc,Pi,Qi,linec,Sflujo,PQ,PV,Slack,Haa,Vcomp,Vmin,Vmax)

    










