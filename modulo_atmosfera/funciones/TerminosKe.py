""" 
Calculo términos de evolución de energía cinética de las perturbaciones
Leandro Diaz, Daniela Risaro
2022

"""

def TermEvolKe(Ke,EstBas_u,EstBas_v,Anomzon_u,Anomzon_v,etae,H):

    """
    Calculo términos de evolución de energía cinética de las perturbaciones
    INPUTS
    Ke: Energía cinética de las perturbaciones    
    EstBas_u: Estado básico de la componente zonal del viento
    EstBas_v: Estado básico de la componente meridional del viento
    Anomzon_u: Anomalía zonal de la componente zonal del viento
    Anomzon_v: Anomalía zonal de la componente meridional del viento   
    etae: Anomalías zonales de la Altura de la superficie libre 
    H: Altura de referencia (convertir el geopotencial en altura)
    OUTPUTS
    advKe: Advección total de energía cinética de las perturbaciones 
    ConvBarot: Conversion barotropica
    DispKe: Dispersión energía cinética
    ConvBaroc: Conversión baroclínica
    """

    #Cargamos librerias necesarias
    import xarray as xr
    import numpy as np

    #Definimos parametros
    Rt=6370*1000 #Radio terrestre [m] 
    #Calculo paso espacial meridional (en metros)
    Per_Grad=np.pi*Rt/180 #Perimetro terrestre por grado de latitud
    dy=Ke['lat'].differentiate('lat')*Per_Grad
    #Calculo paso espacial zonal (en metros). Hay que corregir después por latitud
    dx=Ke['lon'].differentiate('lon')*Per_Grad        
    # Calculo Advección total de energía cinética de las perturbaciones 
    advKe=-(EstBas_u*Ke.differentiate('lon')/(dx*np.cos(Ke['lat']*np.pi/180)) + EstBas_v*Ke.differentiate('lat')/dy + Anomzon_u*Ke.differentiate('lon')/(dx*np.cos(Ke['lat']*np.pi/180)) + Anomzon_v*Ke.differentiate('lat')/dy)
    # Calculo Conversion barotropica
    ConvBarot=-1*H*(np.power(Anomzon_u,2)*EstBas_u.differentiate('lon')/(dx*np.cos(Ke['lat']*np.pi/180))+Anomzon_u*Anomzon_v*EstBas_u['lat'].differentiate('lat')/dy+Anomzon_u*Anomzon_v*EstBas_v.differentiate('lon')/(dx*np.cos(Ke['lat']*np.pi/180))+np.power(Anomzon_v,2)*EstBas_v.differentiate('lat')/dy)
    # Calculo Dispersión energía cinética
    DispKe=-9.8*H*((etae*Anomzon_u).differentiate('lon')/(dx*np.cos(Ke['lat']*np.pi/180))+(etae*Anomzon_v).differentiate('lat')/dy)
    #  Calculo Conversión baroclínica
    ConvBaroc=9.8*H*etae*(Anomzon_u.differentiate('lon')/(dx*np.cos(Ke['lat']*np.pi/180))+Anomzon_v.differentiate('lat')/dy)
    return advKe, ConvBarot, DispKe, ConvBaroc
