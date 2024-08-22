""" 
Cálculo del laplaciano
Leandro Díaz, Daniela Risaro
2022

INPUTS
array_2D: numpy array de dos dimensiones
ds: paso espacial 

"""

def Calc_del2(array_2D, ds):

    #Cargamos las librerias necesarias
    import numpy as np   
    del2 = np.zeros(array_2D.shape, float)
    del2[1:-1, 1:-1] = (array_2D[1:-1,2:] + array_2D[1:-1,:-2] + array_2D[2:,1:-1] + array_2D[:-2,1:-1] - 4.*array_2D[1:-1,1:-1])/(ds*ds)
    return del2