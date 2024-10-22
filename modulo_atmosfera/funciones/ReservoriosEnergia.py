""" 
Calculo Reservorios de Energía
Leandro Diaz, Daniela Risaro
2022

"""

def TermEner(etaz,etae,EstBas_u,EstBas_v,Anomzon_u,Anomzon_v,H):

    """
    Calculo Reservorios de energía
    INPUTS
    etaz: Altura de la superficie libre del estado básico
    etae: Anomalías zonales de la Altura de la superficie libre 
    EstBas_u: Estado básico de la componente zonal del viento
    EstBas_v: Estado básico de la componente meridional del viento
    Anomzon_u: Anomalía zonal de la componente zonal del viento
    Anomzon_v: Anomalía zonal de la componente meridional del viento   
    H: Altura de referencia (convertir el geopotencial en altura)
    OUTPUTS
    Az: Energía potencial del estado básico
    Kz: Energía cinética del estado básico
    Ae: Energía potencial de las perturbaciones
    Ke: Energía cinética de las perturbaciones
    """

    #Cargamos librerias necesarias
    import xarray as xr
    import numpy as np
    
    #Definimos parametros
    rho=1 #Densidad [kg/m3]
    g=9.8 #Aceleración de la gravedad [M/s2]
    #Energia potencial estado basico
    Az = rho*(g/2.)*(np.power(etaz,2)+np.power(H,2))
    #Energia cinetica estado basico
    Kz = rho*(H/2.)*(np.power(EstBas_u,2)+np.power(EstBas_v,2)) 
    #Energia potencial perturbaciones
    Ae = rho*(g/2.)*(np.power(etae,2))
    #Energia cinetica perturbaciones
    Ke = rho*(H/2.)*(np.power(Anomzon_u,2)+np.power(Anomzon_v,2)) 
    return Az, Kz, Ae, Ke
