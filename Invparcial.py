# -*- coding: utf-8 -*-
"""
Created on Sat Aug  7 23:44:50 2021

@author: Orlando RB
"""

# Construcción de la matriz de admitancias
def shipley(B,pos):
    import numpy as np
    n = len(B);        # Longitud
    A = np.zeros((n,n),dtype=complex)    # Se crea una matriz de ceros del tamaño de la entrada
    a_kk = B[pos,pos]      # Variable auxiliar, que me indica el valor de la diagonal que voy a reemplazar
    b_aux = np.zeros(n)   # Vector auxiliar para realizar los calculos de sumas
    # Cambiamos la fila que se desea pivotear por su equivalente al realizar el cambio 
    b_aux = b_aux - B[pos,:]/a_kk
    A[pos,:] = 1*b_aux
    A[pos,pos] = A[pos,pos]/a_kk
    # A(pos,:)=b_aux'; 
    # Calculamos  El resto de sus filas
    for i in range(0,n):
        if i != pos:
            A[i,:] = B[i,:] + (B[i,pos]*b_aux);
            A[i,pos] = -B[i,pos]/a_kk;
    return A