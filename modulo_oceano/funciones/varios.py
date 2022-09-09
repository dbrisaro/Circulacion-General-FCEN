"""
Operaciones para resolver la P1 del modulo de Oceano
Leandro Diaz y Daniela Risaro
2022

"""

import math

# Escalas tipicas
tau = 0.25          #Tension del viento [N/m^2]
L = 4000000         # Longitud de la cuenca [m]
D = 2500            # Profundidad [m]
beta = 2e-11        # Coeficiente de Coriolis [1/(s*m)]
rho = 1025          # Densidad [kg/m^3]


# Parametros para la dimensionalizacion
U = (2*math.pi*tau)/(rho*D*beta*L)                              # Velocidad
Ro = (2*math.pi*tau)/(rho*D*(math.pow(beta,2)*(math.pow(L,3)))) # Numero de Rossby

# Energia cinetica
tke =  QG_diag[:,3]  # la tercera columna tiene la Energia cinetica total 


# Dimensionalizacion de las variables
psidim = psiadim*U*L
vortdim = vortadim*U/L

# Transporte meridional promediado en la vertical (derivada zonal de la función corriente multiplicada por la profundidad)
trans_mer = np.diff(psiadim, n=1, axis=1)*D

# Terminos ecuacion de stommel
ter1 = np.diff(psiadim, n=1, axis=1)
ter1_LatCent = np.squeeze(ter1[int(np.size(ter1,0)/2),:])*Lx/(nx-1)

ter2 = -QG_curlw[int(np.size(ter1,0)/2),1:(nx+1)]
ter3 = Ef*vortadim[int(np.size(ter1,0)/2),:]

# Terminos ecuacion Munk
ds = 0.1
term1 = np.diff(psiadim,n=1,axis=1)[int(np.size(psiadim,0)/2),:]/ds
term2 = -QG_curlw[int(np.size(QG_curlw,0)/2),1:-1]    
term3 = -Ev1*Calc_del2(vortadim,ds)[int(np.size(vortadim,0)/2),:]  

