""" 
Calculo Anomalías zonales de función corriente
Leandro Diaz, Daniela Risaro
2022

"""

def Anomzon_Stream(stream,lon):

    """
    Calculo Anomalías zonales de función corriente
    Calculo Estado básico función corriente
    INPUTS
    lon: longitudes
    stream: Campos de función corriente (todos los tiempos)
    OUTPUTS
    AnomStream: Anomalías zonales de función corriente
    BasicStream: Estado básico función corriente
    """

    #Cargamos librerias necesarias
    import math
    import numpy as np
    
    #Calculamos el Promedio zonal de la función corriente para el estado Básico (previo al restart)
    l,m=np.meshgrid(lon,np.mean(stream[49,:,:],axis=1))
    #Calculamos el Promedio zonal de la función corriente promedio de la simulación después del restart (agregando la nueva perturbación)
    n=np.mean(stream[50:60,:,:],axis=0)
    l,n1=np.meshgrid(lon,np.mean(n,axis=1))
    #Calculamos el estado básico a partir de la correción del estado básico
    BasicStream=stream[49,:,:]-m+n1
    AnomStream=stream[50:60,:,:] #Función corriente después del restart
    i=0
    while i<np.size(AnomStream,0):
        AnomStream[i,:,:]=AnomStream[i,:,:]-basicstream
        i=i+1
    return AnomStream, BasicStream