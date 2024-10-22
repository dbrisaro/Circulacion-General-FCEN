""" 
Calculo promedio global pesado arealmente 
Leandro Diaz, Daniela Risaro
2023

"""

def PromGlobal(var):

    """
    Calculo sumatoria global pesado arealmente 
    INPUTS
    var: Array de variables para la simulación (todos los tiempos)
    OUTPUTS
    PromGlob: Promedio global pesado arealmente 
    """

    #Cargamos librerias necesarias
    import xarray as xr
    import numpy as np

    #Calculamos pesos para hacer el calculo pesado por área
    weights = np.cos(np.deg2rad(var.lat))
    #Asignamos los pesos
    var_weighted = var.weighted(weights)
    #Calculamos el promedio pesado por área
    var_mean = var_weighted.mean(("lon", "lat"))
    
    return var_mean
