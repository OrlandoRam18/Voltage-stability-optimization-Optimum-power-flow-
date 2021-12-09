# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 16:32:16 2020

@author: orlan_qtpj2q1
"""


def UpdateDesitionVariablesA(Vectorobjetivo, vec1, vec2):
    import numpy as np
    indi=Vectorobjetivo.argmin() # index of the best solution at this iteration
    # Best of each component
    Bestxd1=vec1[indi];    Bestxd2=vec2[indi];
    # Find the worst result
    indo=Vectorobjetivo.argmax() # index of the best solution at this iteration
    # Worst of each component
    Worstxd1=vec1[indo]; Worstxd2=vec2[indo]
    # Nes values of the desition variables
    X1=np.zeros((len(vec1),len(np.transpose(vec1)))); X2=np.zeros((len(vec2),len(np.transpose(vec2))))
    for k in range(0,len(np.transpose(vec1))):
        for i in range(0,len(vec1)):
            r=np.random.rand(2);
            X1[i,k]=vec1[i,k]+r[0]*(Bestxd1[k]-abs(vec1[i,k]))-r[1]*(Worstxd1[k]-abs(vec1[i,k]))
    for k in range(0,len(np.transpose(vec2))):
        for i in range(0,len(vec2)):
            r=np.random.rand(2);
            X2[i,k]=vec2[i,k]+r[0]*(Bestxd2[k]-abs(vec2[i,k]))-r[1]*(Worstxd2[k]-abs(vec2[i,k]))
    return X1, X2



def UpdateDesitionVariablesB(Vectorobjetivo, vec1, vec2, vec3):
    import numpy as np
    indi=Vectorobjetivo.argmin() # index of the best solution at this iteration
    # Best of each component
    Bestxd1=vec1[indi];    Bestxd2=vec2[indi]; Bestxd3=vec3[indi]
    # Find the worst result
    indo=Vectorobjetivo.argmax() # index of the best solution at this iteration
    # Worst of each component
    Worstxd1=vec1[indo]; Worstxd2=vec2[indo]; Worstxd3=vec3[indo]
    # Nes values of the desition variables
    X1=np.zeros((len(vec1),len(np.transpose(vec1)))); X2=np.zeros((len(vec2),len(np.transpose(vec2))))
    X3=np.zeros((len(vec3),len(np.transpose(vec3))))
    for k in range(0,len(np.transpose(vec1))):
        for i in range(0,len(vec1)):
            r=np.random.rand(2);
            X1[i,k]=vec1[i,k]+r[0]*(Bestxd1[k]-abs(vec1[i,k]))-r[1]*(Worstxd1[k]-abs(vec1[i,k]))
    for k in range(0,len(np.transpose(vec2))):
        for i in range(0,len(vec2)):
            r=np.random.rand(2);
            X2[i,k]=vec2[i,k]+r[0]*(Bestxd2[k]-abs(vec2[i,k]))-r[1]*(Worstxd2[k]-abs(vec2[i,k]))
    for k in range(0,len(np.transpose(vec3))):
        for i in range(0,len(vec3)):
            r=np.random.rand(2);
            X3[i,k]=vec3[i,k]+r[0]*(Bestxd3[k]-abs(vec3[i,k]))-r[1]*(Worstxd3[k]-abs(vec3[i,k]))

    return X1, X2, X3



def UpdateDesitionVariablesC(Vectorobjetivo, vec1, vec2, vec3, vec4, vec5, vec6):
    import numpy as np
    indi=Vectorobjetivo.argmin() # index of the best solution at this iteration
    # Best of each component
    Bestxd1=vec1[indi];    Bestxd2=vec2[indi]; Bestxd3=vec3[indi]; Bestxd4=vec4[indi]
    Bestxd5=vec5[indi]; Bestxd6=vec6[indi]
    # Find the worst result
    indo=Vectorobjetivo.argmax() # index of the best solution at this iteration
    # Worst of each component
    Worstxd1=vec1[indo]; Worstxd2=vec2[indo]; Worstxd3=vec3[indo]; Worstxd4=vec4[indo]
    Worstxd5=vec5[indo]; Worstxd6=vec6[indo]
    # Nes values of the desition variables
    X1=np.zeros((len(vec1),len(np.transpose(vec1)))); X2=np.zeros((len(vec2),len(np.transpose(vec2))))
    X3=np.zeros((len(vec3),len(np.transpose(vec3)))); X4=np.zeros((len(vec4),len(np.transpose(vec4))))
    X5=np.zeros((len(vec5),len(np.transpose(vec5)))); X6=np.zeros((len(vec6),len(np.transpose(vec6))))
    for k in range(0,len(np.transpose(vec1))):
        for i in range(0,len(vec1)):
            r=np.random.rand(2);
            X1[i,k]=vec1[i,k]+r[0]*(Bestxd1[k]-abs(vec1[i,k]))-r[1]*(Worstxd1[k]-abs(vec1[i,k]))
    for k in range(0,len(np.transpose(vec2))):
        for i in range(0,len(vec2)):
            r=np.random.rand(2);
            X2[i,k]=vec2[i,k]+r[0]*(Bestxd2[k]-abs(vec2[i,k]))-r[1]*(Worstxd2[k]-abs(vec2[i,k]))
    for k in range(0,len(np.transpose(vec3))):
        for i in range(0,len(vec3)):
            r=np.random.rand(2);
            X3[i,k]=vec3[i,k]+r[0]*(Bestxd3[k]-abs(vec3[i,k]))-r[1]*(Worstxd3[k]-abs(vec3[i,k]))
    for k in range(0,len(np.transpose(vec4))):
        for i in range(0,len(vec4)):
            r=np.random.rand(2);
            X4[i,k]=vec4[i,k]+r[0]*(Bestxd4[k]-abs(vec4[i,k]))-r[1]*(Worstxd4[k]-abs(vec4[i,k]))
    for k in range(0,len(np.transpose(vec5))):
        for i in range(0,len(vec5)):
            r=np.random.rand(2);
            X5[i,k]=vec5[i,k]+r[0]*(Bestxd5[k]-abs(vec5[i,k]))-r[1]*(Worstxd5[k]-abs(vec5[i,k]))
    for k in range(0,len(np.transpose(vec6))):
        for i in range(0,len(vec6)):
            r=np.random.rand(2);
            X6[i,k]=vec6[i,k]+r[0]*(Bestxd6[k]-abs(vec6[i,k]))-r[1]*(Worstxd6[k]-abs(vec6[i,k]))
    return X1, X2, X3, X4, X5, X6



def UpdateDesitionVariablesD(Vectorobjetivo, vec1, vec2, vec3, vec4, vec5):
    import numpy as np
    indi=Vectorobjetivo.argmin() # index of the best solution at this iteration
    # Best of each component
    Bestxd1=vec1[indi];    Bestxd2=vec2[indi]; Bestxd3=vec3[indi]; Bestxd4=vec4[indi]
    Bestxd5=vec5[indi]
    # Find the worst result
    indo=Vectorobjetivo.argmax() # index of the best solution at this iteration
    # Worst of each component
    Worstxd1=vec1[indo]; Worstxd2=vec2[indo]; Worstxd3=vec3[indo]; Worstxd4=vec4[indo]
    Worstxd5=vec5[indo]
    # Nes values of the desition variables
    X1=np.zeros((len(vec1),len(np.transpose(vec1)))); X2=np.zeros((len(vec2),len(np.transpose(vec2))))
    X3=np.zeros((len(vec3),len(np.transpose(vec3)))); X4=np.zeros((len(vec4),len(np.transpose(vec4))))
    X5=np.zeros((len(vec5),len(np.transpose(vec5))))
    for k in range(0,len(np.transpose(vec1))):
        for i in range(0,len(vec1)):
            r=np.random.rand(2);
            X1[i,k]=vec1[i,k]+r[0]*(Bestxd1[k]-abs(vec1[i,k]))-r[1]*(Worstxd1[k]-abs(vec1[i,k]))
    for k in range(0,len(np.transpose(vec2))):
        for i in range(0,len(vec2)):
            r=np.random.rand(2);
            X2[i,k]=vec2[i,k]+r[0]*(Bestxd2[k]-abs(vec2[i,k]))-r[1]*(Worstxd2[k]-abs(vec2[i,k]))
    for k in range(0,len(np.transpose(vec3))):
        for i in range(0,len(vec3)):
            r=np.random.rand(2);
            X3[i,k]=vec3[i,k]+r[0]*(Bestxd3[k]-abs(vec3[i,k]))-r[1]*(Worstxd3[k]-abs(vec3[i,k]))
    for k in range(0,len(np.transpose(vec4))):
        for i in range(0,len(vec4)):
            r=np.random.rand(2);
            X4[i,k]=vec4[i,k]+r[0]*(Bestxd4[k]-abs(vec4[i,k]))-r[1]*(Worstxd4[k]-abs(vec4[i,k]))
    for k in range(0,len(np.transpose(vec5))):
        for i in range(0,len(vec5)):
            r=np.random.rand(2);
            X5[i,k]=vec5[i,k]+r[0]*(Bestxd5[k]-abs(vec5[i,k]))-r[1]*(Worstxd5[k]-abs(vec5[i,k]))

    return X1, X2, X3, X4, X5


def UpdateDesitionVariablesOmega(Vectorobjetivo, vec1, vec2, vec3, vec4):
    import numpy as np
    indi=Vectorobjetivo.argmin() # index of the best solution at this iteration
    # Best of each component
    Bestxd1=vec1[indi];    Bestxd2=vec2[indi]; Bestxd3=vec3[indi]; Bestxd4=vec4[indi]
    # Find the worst result
    indo=Vectorobjetivo.argmax() # index of the best solution at this iteration
    # Worst of each component
    Worstxd1=vec1[indo]; Worstxd2=vec2[indo]; Worstxd3=vec3[indo]; Worstxd4=vec4[indo]
    # Nes values of the desition variables
    X1=np.zeros((len(vec1),len(np.transpose(vec1)))); X2=np.zeros((len(vec2),len(np.transpose(vec2))))
    X3=np.zeros((len(vec3),len(np.transpose(vec3)))); X4=np.zeros((len(vec4),len(np.transpose(vec4))))
    for k in range(0,len(np.transpose(vec1))):
        for i in range(0,len(vec1)):
            r=np.random.rand(2);
            X1[i,k]=vec1[i,k]+r[0]*(Bestxd1[k]-abs(vec1[i,k]))-r[1]*(Worstxd1[k]-abs(vec1[i,k]))
    for k in range(0,len(np.transpose(vec2))):
        for i in range(0,len(vec2)):
            r=np.random.rand(2);
            X2[i,k]=vec2[i,k]+r[0]*(Bestxd2[k]-abs(vec2[i,k]))-r[1]*(Worstxd2[k]-abs(vec2[i,k]))
    for k in range(0,len(np.transpose(vec3))):
        for i in range(0,len(vec3)):
            r=np.random.rand(2);
            X3[i,k]=vec3[i,k]+r[0]*(Bestxd3[k]-abs(vec3[i,k]))-r[1]*(Worstxd3[k]-abs(vec3[i,k]))
    for k in range(0,len(np.transpose(vec4))):
        for i in range(0,len(vec4)):
            r=np.random.rand(2);
            X4[i,k]=vec4[i,k]+r[0]*(Bestxd4[k]-abs(vec4[i,k]))-r[1]*(Worstxd4[k]-abs(vec4[i,k]))
    return X1, X2, X3, X4


def UpdateDesitionVariablesAA(Vectorobjetivo, vec1, vec2):
    import numpy as np
    indi=Vectorobjetivo.argmin() # index of the best solution at this iteration
    # Best of each component
    Bestxd1=vec1[indi];    Bestxd2=vec2[indi];
    # Find the worst result
    indo=Vectorobjetivo.argmax() # index of the best solution at this iteration
    # Worst of each component
    Worstxd1=vec1[indo]; Worstxd2=vec2[indo]
    # Nes values of the desition variables
    X1=np.zeros((len(vec1),len(np.transpose(vec1)))); X2=np.zeros((len(vec2),len(np.transpose(vec2))))
    for k in range(0,len(np.transpose(vec1))):
        for i in range(0,len(vec1)):
            r=np.random.rand(2);
            if r[0] < 0.7: r[1] =r[0]/0.7
            elif r[0] >=0.7: r[1] = (10/3)*(1-r[0])
            X1[i,k]=vec1[i,k]+r[0]*(Bestxd1[k]-abs(vec1[i,k]))-r[1]*(Worstxd1[k]-abs(vec1[i,k]))
    for k in range(0,len(np.transpose(vec2))):
        for i in range(0,len(vec2)):
            if r[0] < 0.7: r[1] =r[0]/0.7
            elif r[0] >=0.7: r[1] = (10/3)*(1-r[0])
            r=np.random.rand(2);
            X2[i,k]=vec2[i,k]+r[0]*(Bestxd2[k]-abs(vec2[i,k]))-r[1]*(Worstxd2[k]-abs(vec2[i,k]))
    return X1, X2



def UpdateDesitionVariablesBB(Vectorobjetivo, vec1, vec2, vec3):
    import numpy as np
    indi=Vectorobjetivo.argmin() # index of the best solution at this iteration
    # Best of each component
    Bestxd1=vec1[indi];    Bestxd2=vec2[indi]; Bestxd3=vec3[indi]
    # Find the worst result
    indo=Vectorobjetivo.argmax() # index of the best solution at this iteration
    # Worst of each component
    Worstxd1=vec1[indo]; Worstxd2=vec2[indo]; Worstxd3=vec3[indo]
    # Nes values of the desition variables
    X1=np.zeros((len(vec1),len(np.transpose(vec1)))); X2=np.zeros((len(vec2),len(np.transpose(vec2))))
    X3=np.zeros((len(vec3),len(np.transpose(vec3))))
    for k in range(0,len(np.transpose(vec1))):
        for i in range(0,len(vec1)):
            r=np.random.rand(2);
            if r[0] < 0.7: r[1] =r[0]/0.7
            elif r[0] >=0.7: r[1] = (10/3)*(1-r[0])
            X1[i,k]=vec1[i,k]+r[0]*(Bestxd1[k]-abs(vec1[i,k]))-r[1]*(Worstxd1[k]-abs(vec1[i,k]))
    for k in range(0,len(np.transpose(vec2))):
        for i in range(0,len(vec2)):
            r=np.random.rand(2);
            if r[0] < 0.7: r[1] =r[0]/0.7
            elif r[0] >=0.7: r[1] = (10/3)*(1-r[0])
            X2[i,k]=vec2[i,k]+r[0]*(Bestxd2[k]-abs(vec2[i,k]))-r[1]*(Worstxd2[k]-abs(vec2[i,k]))
    for k in range(0,len(np.transpose(vec3))):
        for i in range(0,len(vec3)):
            r=np.random.rand(2);
            if r[0] < 0.7: r[1] =r[0]/0.7
            elif r[0] >=0.7: r[1] = (10/3)*(1-r[0])
            X3[i,k]=vec3[i,k]+r[0]*(Bestxd3[k]-abs(vec3[i,k]))-r[1]*(Worstxd3[k]-abs(vec3[i,k]))
    return X1, X2, X3



def UpdateDesitionVariablesCC(Vectorobjetivo, vec1, vec2, vec3, vec4):
    import numpy as np
    indi=Vectorobjetivo.argmin() # index of the best solution at this iteration
    # Best of each component
    Bestxd1=vec1[indi];    Bestxd2=vec2[indi]; Bestxd3=vec3[indi]; Bestxd4=vec4[indi]
    # Find the worst result
    indo=Vectorobjetivo.argmax() # index of the best solution at this iteration
    # Worst of each component
    Worstxd1=vec1[indo]; Worstxd2=vec2[indo]; Worstxd3=vec3[indo]; Worstxd4=vec4[indo]
    # Nes values of the desition variables
    X1=np.zeros((len(vec1),len(np.transpose(vec1)))); X2=np.zeros((len(vec2),len(np.transpose(vec2))))
    X3=np.zeros((len(vec3),len(np.transpose(vec3)))); X4=np.zeros((len(vec4),len(np.transpose(vec4))))
    for k in range(0,len(np.transpose(vec1))):
        for i in range(0,len(vec1)):
            r=np.random.rand(2);
            if r[0] < 0.7: r[1] =r[0]/0.7
            elif r[0] >=0.7: r[1] = (10/3)*(1-r[0])
            X1[i,k]=vec1[i,k]+r[0]*(Bestxd1[k]-abs(vec1[i,k]))-r[1]*(Worstxd1[k]-abs(vec1[i,k]))
    for k in range(0,len(np.transpose(vec2))):
        for i in range(0,len(vec2)):
            r=np.random.rand(2);
            if r[0] < 0.7: r[1] =r[0]/0.7
            elif r[0] >=0.7: r[1] = (10/3)*(1-r[0])
            X2[i,k]=vec2[i,k]+r[0]*(Bestxd2[k]-abs(vec2[i,k]))-r[1]*(Worstxd2[k]-abs(vec2[i,k]))
    for k in range(0,len(np.transpose(vec3))):
        for i in range(0,len(vec3)):
            r=np.random.rand(2);
            if r[0] < 0.7: r[1] =r[0]/0.7
            elif r[0] >=0.7: r[1] = (10/3)*(1-r[0])
            X3[i,k]=vec3[i,k]+r[0]*(Bestxd3[k]-abs(vec3[i,k]))-r[1]*(Worstxd3[k]-abs(vec3[i,k]))
    for k in range(0,len(np.transpose(vec4))):
        for i in range(0,len(vec4)):
            r=np.random.rand(2);
            if r[0] < 0.7: r[1] =r[0]/0.7
            elif r[0] >=0.7: r[1] = (10/3)*(1-r[0])
            X4[i,k]=vec4[i,k]+r[0]*(Bestxd4[k]-abs(vec4[i,k]))-r[1]*(Worstxd4[k]-abs(vec4[i,k]))
    return X1, X2, X3, X4



def UpdateDesitionVariablesDD(Vectorobjetivo, vec1, vec2, vec3, vec4, vec5):
    import numpy as np
    indi=Vectorobjetivo.argmin() # index of the best solution at this iteration
    # Best of each component
    Bestxd1=vec1[indi];    Bestxd2=vec2[indi]; Bestxd3=vec3[indi]; Bestxd4=vec4[indi]
    Bestxd5=vec5[indi]
    # Find the worst result
    indo=Vectorobjetivo.argmax() # index of the best solution at this iteration
    # Worst of each component
    Worstxd1=vec1[indo]; Worstxd2=vec2[indo]; Worstxd3=vec3[indo]; Worstxd4=vec4[indo]
    Worstxd5=vec5[indo]
    # Nes values of the desition variables
    X1=np.zeros((len(vec1),len(np.transpose(vec1)))); X2=np.zeros((len(vec2),len(np.transpose(vec2))))
    X3=np.zeros((len(vec3),len(np.transpose(vec3)))); X4=np.zeros((len(vec4),len(np.transpose(vec4))))
    X5=np.zeros((len(vec5),len(np.transpose(vec5))))
    for k in range(0,len(np.transpose(vec1))):
        for i in range(0,len(vec1)):
            r=np.random.rand(2);
            if r[0] < 0.7: r[1] =r[0]/0.7
            elif r[0] >=0.7: r[1] = (10/3)*(1-r[0])
            X1[i,k]=vec1[i,k]+r[0]*(Bestxd1[k]-abs(vec1[i,k]))-r[1]*(Worstxd1[k]-abs(vec1[i,k]))
    for k in range(0,len(np.transpose(vec2))):
        for i in range(0,len(vec2)):
            r=np.random.rand(2);
            if r[0] < 0.7: r[1] =r[0]/0.7
            elif r[0] >=0.7: r[1] = (10/3)*(1-r[0])
            X2[i,k]=vec2[i,k]+r[0]*(Bestxd2[k]-abs(vec2[i,k]))-r[1]*(Worstxd2[k]-abs(vec2[i,k]))
    for k in range(0,len(np.transpose(vec3))):
        for i in range(0,len(vec3)):
            r=np.random.rand(2);
            if r[0] < 0.7: r[1] =r[0]/0.7
            elif r[0] >=0.7: r[1] = (10/3)*(1-r[0])
            X3[i,k]=vec3[i,k]+r[0]*(Bestxd3[k]-abs(vec3[i,k]))-r[1]*(Worstxd3[k]-abs(vec3[i,k]))
    for k in range(0,len(np.transpose(vec4))):
        for i in range(0,len(vec4)):
            r=np.random.rand(2);
            if r[0] < 0.7: r[1] =r[0]/0.7
            elif r[0] >=0.7: r[1] = (10/3)*(1-r[0])
            X4[i,k]=vec4[i,k]+r[0]*(Bestxd4[k]-abs(vec4[i,k]))-r[1]*(Worstxd4[k]-abs(vec4[i,k]))
    for k in range(0,len(np.transpose(vec5))):
        for i in range(0,len(vec5)):
            r=np.random.rand(2);
            if r[0] < 0.7: r[1] =r[0]/0.7
            elif r[0] >=0.7: r[1] = (10/3)*(1-r[0])
            X5[i,k]=vec5[i,k]+r[0]*(Bestxd5[k]-abs(vec5[i,k]))-r[1]*(Worstxd5[k]-abs(vec5[i,k]))
    return X1, X2, X3, X4, X5


def UpdateDesitionVariables(Vectorobjetivo, vec1, vec2, vec3):
    import numpy as np
    indi=Vectorobjetivo.argmin() # index of the best solution at this iteration
    # Best of each component
    Bestxd1=vec1[indi];    Bestxd2=vec2[indi]; Bestxd3=vec3[indi]
    # Find the worst result
    indo=Vectorobjetivo.argmax() # index of the best solution at this iteration
    # Worst of each component
    Worstxd1=vec1[indo]; Worstxd2=vec2[indo]; Worstxd3=vec3[indo]
    # Nes values of the desition variables
    X1=np.zeros((len(vec1),len(np.transpose(vec1)))); X2=np.zeros((len(vec2),len(np.transpose(vec2))))
    X3=np.zeros((len(vec3),len(np.transpose(vec3))))
    for k in range(0,len(np.transpose(vec1))):
        for i in range(0,len(vec1)):
            r=np.random.rand(2);
            X1[i,k]=vec1[i,k]+r[0]*(Bestxd1[k]-abs(vec1[i,k]))-r[1]*(Worstxd1[k]-abs(vec1[i,k]))
    for k in range(0,len(np.transpose(vec2))):
        for i in range(0,len(vec2)):
            r=np.random.rand(2);
            X2[i,k]=vec2[i,k]+r[0]*(Bestxd2[k]-abs(vec2[i,k]))-r[1]*(Worstxd2[k]-abs(vec2[i,k]))
    for k in range(0,len(np.transpose(vec3))):
        for i in range(0,len(vec3)):
            r=np.random.rand(2);
            X3[i,k]=vec3[i,k]+r[0]*(Bestxd3[k]-abs(vec3[i,k]))-r[1]*(Worstxd3[k]-abs(vec3[i,k]))
    return X1, X2, X3