""" 
Calculo Anomalias zonales de una variable según corección estado básico
Leandro Diaz, Daniela Risaro
2023

"""

def Anomzon(var):

    """
    Calculo Anomalías zonales
    Calculo Estado básico
    INPUTS
    var: Array de variables para la simulación (todos los tiempos)
    OUTPUTS
    AnomZonVar: Anomalías zonales
    EstBasVar: Estado básico
    """

    #Cargamos librerias necesarias
    import xarray as xr
    
    #Promedio zonal para el estado Básico
    PromZonalEstBas=var.sel(time=50).mean(dim='lon')
    #Promedio zonal para la simulación después del restart
    PromZonalRestart=var.sel(time=slice(51,60)).mean(dim=('lon','time'))
    #Correción con la cual vamos a calcular las anomalías
    #Calculamos estado básico
    EstBasVar=var.sel(time=50)-PromZonalEstBas+PromZonalRestart
    #Calculamos anomalías
    AnomZonVar=var.sel(time=slice(51,60))-EstBasVar
    return AnomZonVar,EstBasVar
