""" 
Calculo Numero de Onda Estacionario (Ks)
Leandro Diaz, Daniela Risaro
2022

"""

def Ks(lat,u,GradMerVortAbs):

    """
    Calculo Numero de Onda Estacionario 
    Calculo Numero de Onda Estacionario Planetario
    INPUTS
    lat: latitudes
    u: Campo de velocidad zonal para el estado estacionario    
    GradMerVortAbs: Gradiente meridional de vorticidad absoluta
    OUTPUTS
    Ks: Numero de Onda Estacionario
    Ks_Plan: Numero de Onda Estacionario Planetario 
    """

    #Cargamos librerias necesarias
    import xarray as xr
    import math
    import numpy as np
    
    #Definimos parametros
    Rt=6370*1000 #Radio terrestre [m]    
    #Calculo Numero de Onda Estacionario
    Ks=np.sqrt(GradMerVortAbs/u)
    #Calculo Numero de Onda Estacionario Planetario
    KsPlan=Rt*np.cos(math.pi*lat/180)*Ks    
    
    return Ks,KsPlan