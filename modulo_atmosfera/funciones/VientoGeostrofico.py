""" 
Calculo Viento Geostrófico
Leandro Diaz, Daniela Risaro
2022

"""

def Ugeo(h):

    """
    Calculo Viento Geostrófico
    INPUTS
    h: Geopotencial de la capa 
    OUTPUTS
    ug: Componente zonal del Viento geostrófico
    vg: Componente meridional del Viento geostrófico
    """

    #Cargamos librerias necesarias
    import xarray as xr
    import numpy as np
    
    #Definimos parametros
    omega=7.29*0.00001 #Frecuencia angular terrestre [1/s]
    Rt=6370*1000 #Radio terrestre [m] 
    #Calculamos el parámetro de Coriolis
    f=2*omega*np.sin(h['lat']*np.pi/180)
    #Calculo paso espacial meridional (en metros)
    Per_Grad=np.pi*Rt/180 #Perimetro terrestre por grado de latitud
    dy=h['lat'].differentiate('lat')*Per_Grad
    #Calculo paso espacial zonal (en metros). Hay que corregir después por latitud
    dx=h['lon'].differentiate('lon')*Per_Grad    
    #Calculo viento geostrófico
    Inv_f=1/f
    ug=-Inv_f*h.differentiate('lat')/dy
    vg=Inv_f*h.differentiate('lon')/(dx*np.cos(h['lat']*np.pi/180))

    return ug, vg
