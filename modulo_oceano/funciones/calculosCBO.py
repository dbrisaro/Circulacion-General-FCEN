""" 
Cálculos CBO
Leandro Díaz, Daniela Risaro
2022

"""

def Calc_TrasMer_CBO_LatCent(X,psiF,U,L,D):

    """
    Calculo de Extensión de la Corriente de Borde Oeste en la latitud central,
    Calculo de Transporte total meridional de la Corriente de Borde Oeste en la latitud central
    Calculo de transporte total meridional en la latitud central
    
    INPUTS

    OUTPUTS
    """

    TrasMer=np.diff(psiF,n=1,axis=1)*U*L*D
    TrasMer_LatCent=TrasMer[int(np.size(TrasMer,0)/2),:]
    X_mod=X[0:-1]-(X[1]-X[0])/2

    a=TrasMer_LatCent[0]/abs(TrasMer_LatCent[0])
    i=1
    m=0
    while i<len(TrasMer_LatCent):
        b=TrasMer_LatCent[i]/abs(TrasMer_LatCent[i])
        if b==a:
            i=i+1
        else:
            m=i
            break    
    Limite_CBO_LatCent=X_mod[m]+(X[1]-X[0])/2
    TrasMer_CBO_LatCent=sum(TrasMer_LatCent[0:m+1])/1000000
    TrasMer_total_LatCent=sum(TrasMer_LatCent)/1000000 
    
    return Limite_CBO_LatCent,TrasMer_CBO_LatCent,TrasMer_total_LatCent