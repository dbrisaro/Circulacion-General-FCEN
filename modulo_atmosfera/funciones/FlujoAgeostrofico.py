""" 
Calculo flujo ageostrofico de geopotencial
Leandro Diaz, Daniela Risaro
2022

"""

def FlujoAgeostrof(etae,Anomzon_uag,Anomzon_vag):

    """
    Calculo términos de evolución de energía cinética de las perturbaciones
    INPUTS
    etae: Anomalías zonales de la Altura de la superficie libre 
    Anomzon_uag: Anomalía zonal de la componente zonal del viento ageostrofico
    Anomzon_vag: Anomalía zonal de la componente meridional del viento ageostrofico
    OUTPUTS
    Fagx: Componente zonal del flujo ageostrófico de geopotencial
    Fagy: Componente meridional del flujo ageostrófico de geopotencial
    """

    #Cargamos librerias necesarias
    import xarray as xr
    
    #Calculamos el flujo ageostrofico de geopotencial
    Fagx=etae*Anomzon_uag
    Fagy=etae*Anomzon_vag
    
    return Fagx, Fagy
