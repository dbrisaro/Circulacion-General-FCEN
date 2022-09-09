""" 
Número de iteraciones para alcanzar el estado estacionario
Leandro Díaz, Daniela Risaro
2022

"""
def calc_tiempo_estabilizacion(tke):
    ECfinal=tke[-1]
    error_ECfinal=ECfinal/100         
    i=1
    while i<len(tke):    
        if abs(tke[-i]-ECfinal)<error_ECfinal:
            i=i+1
        else:
            TiempoEst=len(tke)-i
            break  
    return TiempoEst