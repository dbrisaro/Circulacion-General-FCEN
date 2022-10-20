""" 
Calculo Gradientes meridionales de vorticidad relativa, planetaria y absoluta considerando esfericidad terrerstre
Leandro Diaz, Daniela Risaro
2022

"""

def Grad_Med_Vort(lat,u):

    """
    Calculo Gradiente meridional de vorticidad relativa
    Calculo Gradiente meridional de vorticidad planetaria (beta)
    Calculo Gradiente meridional de vorticidad absoluta 
    INPUTS
    lat: latitudes
    u: Campo de velocidad zonal para el estado estacionario
    OUTPUTS
    GradMerVortRel: Gradiente meridional de vorticidad relativa
    GradMerVortPlan: Gradiente meridional de vorticidad planetaria (beta)
    GradMerVortAbs: Gradiente meridional de vorticidad absoluta    
    """

    #Cargamos librerias necesarias
    import xarray as xr
    import math
    import numpy as np
    
    #Definimos parametros
    omega=7.29*0.00001 #Frecuencia angular terrestre [1/s]
    Rt=6370*1000 #Radio terrestre [m]    
    #Calculo paso espacial meridional (en metros)
    Per_Grad=np.pi*Rt/180 #Perimetro terrestre por grado de latitud
    dy=lat.differentiate('lat')*Per_Grad
    #Calculo gradiente meridional de vorticidad relativa
    GradMerVortRel=-(u.differentiate('lat')/dy).differentiate('lat')/dy
    #Calculo gradiente meridional de vorticidad planetaria (beta)
    GradMerVortPlan=2*omega*np.cos(math.pi*lat/180)/Rt
    #Calculo gradiente meridional de vorticidad absoluta
    GradMerVortAbs=GradMerVortRel+GradMerVortPlan
    
    return GradMerVortRel,GradMerVortPlan,GradMerVortAbs