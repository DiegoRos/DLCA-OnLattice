"""
Código: Regresión Linear
Creador: Diego Rosenberg

Resumen: Calcula los valores de la pendiente y la ordenada al origen de dos listas.
"""

#=======================================================================================
# Funciones
#=======================================================================================
#Función encargada de cacular la pendiente y ordenada al origen de dos listas.
def linearReg(xlist, ylist):
    #Variables base
    n = len(xlist)
    xsum = sum(xlist)
    ysum = sum(ylist)
    xsquaredsum = sum([x**2 for x in xlist])

    xysum =0
    for i in range(len(xlist)): #Si la lista x es mas larga que la lista y habrá un error.
        xysum = (xlist[i]*ylist[i]) + xysum

    #Fórmulas utilizadas para calcular los coeficientes.
    m = ((n*xysum) - (xsum*ysum)) / ((n*xsquaredsum) - (xsum*xsum))
    b = (ysum - (m*xsum))/n

    return m, b

#Este if esta colocado para que lo que está adentro de este if se corra solo si se corre este código directamente.
if __name__ == '__main__':
    print(LinearReg([1,2],[1,2]))