import sympy as sp
import numpy as np
from scipy.signal import argrelmin, argrelmax
import copy as cp
import pandas as pd
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Metodos auxiliares

#CSV
def leerSenial():
    datos = pd.read_csv(__file__ + '\..\cardio_03.csv', header = 0, names=['Tiempo','Senial', 'r'])
    tArr = np.array(datos['Tiempo']) * 60 / 100
    xtRuidoArr = np.array(datos['Senial'])

    return tArr, xtRuidoArr

# Encuentra valor mas cercano al dato en el arreglo
def indiceDelMasCercano(t0, arr):   
    absDifs = cp.deepcopy(arr)
    absDifs = abs(absDifs - t0)
    return np.argmin(absDifs)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Variables globales utiles

tArr, xtRuidoArr = leerSenial()

extInts = np.array([0., 7.5, 18., 19.5, 23.5, 35., 45., 60.])   # Extremos de los intervalos

for i in range(len(extInts)):
    extInts[i] = tArr[indiceDelMasCercano(extInts[i], tArr)]

cantInts = len(extInts) - 1

cantCoefsPol = 5                             # Cantidad de coeficientes de los polinomios

cantIncognitas = cantInts * cantCoefsPol

coefs = []                                  # Todos los coeficientes de los polinomios
for i in range(cantIncognitas):
    coefs.append(sp.Symbol('a' + str(i)))
    
ecsContYDeriv = []                          # Ecuaciones de continuidad y derivabilidad de los polinomios

coefsIncognitas = set()                     # Coeficientes 'libres' de ecsContYDeriv

coefsDespejados = []                        # Coeficientes en funcion de otros

t = sp.Symbol('t')                          # Variable independiente

xLambdaArr = []                             # Lista con las 'subfunciones' de la funcion a trozos


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Inciso 1

# Arma un sistema de ecuaciones para las condiciones de continuidad y derivabilidad de la funcion a trozos y lo resuelve simbolicamente
def resolverEcsContYDeriv():
    
    global ecsContYDeriv
    global coefsDespejados
    global coefsIncognitas

    # Valores de la funcion y su derivada en los extremos izquierdo y derecho:
    xIzq  = 0.
    dxIzq = 0.
    xDer  = 0.
    dxDer = 0.

    # Armar ecuaciones:

    # Continuidad en el extremo izquierdo
    ecCont = sp.Eq(xi(0).subs(t, extInts[0]), xIzq)
    ecsContYDeriv.append(ecCont)

   # Suavidad en el extremo izquierdo
    ecDeriv = sp.Eq(xi(0).diff(t).subs(t, extInts[0]), dxIzq)
    ecsContYDeriv.append(ecDeriv)
    
    for i in range(cantInts - 1):
        # Continuidad entre polinomios
        ecCont = sp.Eq(xi(i).subs(t, extInts[i + 1]), xi(i + 1).subs(t, extInts[i + 1]))
        ecsContYDeriv.append(ecCont)
        
        # Suavidad entre polinomios
        ecDeriv = sp.Eq(xi(i).diff(t).subs(t, extInts[i + 1]), xi(i + 1).diff(t).subs(t, extInts[i + 1]))
        ecsContYDeriv.append(ecDeriv)
        
    # Continuidad en el extremo derecho
    ecCont = sp.Eq(xi(cantInts - 1).subs(t, extInts[cantInts]), xDer)
    ecsContYDeriv.append(ecCont)

    # Suavidad en el extremo derecho
    ecDeriv = sp.Eq(xi(cantInts - 1).diff(t).subs(t, extInts[cantInts]), dxDer)
    ecsContYDeriv.append(ecDeriv)

    # Resolver:
    ecsContYDeriv = sp.solve(ecsContYDeriv, coefs)

    # Coeficientes: 
    coefsDespejados = list(ecsContYDeriv.keys())

    for cd in coefsDespejados:
        coefsIncognitas.update(ecsContYDeriv[cd].free_symbols)
    coefsIncognitas = list(coefsIncognitas)

# Retorna la expresion simbolica del polinomio del intervalo int. Los coeficientes pueden estar en funcion de otros o ser escalares 
def xi(int):
    
    global t
    x = t * 0

    if int != -1:
        for j in range(cantCoefsPol):
            try:
                c = ecsContYDeriv[coefs[cantCoefsPol * int + j]]
            except:
                c = coefs[cantCoefsPol * int + j]
        
            x = c * t ** j + x

    return x

def obtenerIndiceIntervalo(t0):
    
    ind = cantInts - 1
    if not (t0 < extInts[0] or t0 >= extInts[cantInts]):
        ind = 0
        piso = extInts[ind]
        techo = extInts[ind + 1]
        while not (piso <= t0 and t0 < techo):
            ind += 1
            piso = extInts[ind]
            techo = extInts[ind + 1]

    return ind


def minimosCuadrados(tArr, xtRuidoArr):
    
    # Expresion simbolica del gamma en un intervalo
    def gamma(int):
        
        global t

        subsRes = xi(int)
        coefsPresentes = subsRes.free_symbols.difference({t})
        for c in coefsPresentes:
            subsRes = subsRes.subs(c, 0)

        return subsRes

    # Retorna el lambda del phi asociado a un coeficiente en un intervalo
    def phiLambda(coef, int):
        
        global t

        subsRes = xi(int) - gamma(int)
        coefsPresentes = subsRes.free_symbols.difference({t})
        for c in coefsPresentes:
            if coef == c:
                subsRes = subsRes.subs(coef, 1)
            else:
                subsRes = subsRes.subs(c, 0)

        return sp.lambdify(t, subsRes)


    # Matriz de expresiones lambda de los phis para cada intervalo:
    m = len(coefsIncognitas)
    lambdasPhis = []
    for i in range(m):
        phiInts = []
        for j in range(cantInts):
            phiInts.append(phiLambda(coefsIncognitas[i], j))
        lambdasPhis.append(phiInts)

    # Arreglo con expresiones lambda del gamma para cada intervalo:
    lambdasGamma = []
    for j in range(cantInts):
        lambdasGamma.append(sp.lambdify(t, gamma(j)))
    
    # Crear e inicializar matriz con los valores de los phis en todos los puntos de la señal:
    matrizPhis = np.zeros((m, len(tArr)))

    for i in range(m):
        for j in range(len(tArr)):
            matrizPhis[i][j] = lambdasPhis[i][obtenerIndiceIntervalo(tArr[j])](tArr[j])

    # Crear e inicializar arreglo con los valores del gamma en todos los puntos de la señal:
    gammaArr = np.zeros(len(tArr))

    for j in range(len(tArr)):
        gammaArr[j] = lambdasGamma[obtenerIndiceIntervalo(tArr[j])](tArr[j])

    # Crear e inicializar matriz y termino independiente para minimos cuadrados: 
    matrizMinCuad = np.zeros((m, m))
    terIndMinCuad = np.zeros(m)

    for i in range(m):
        for j in range(i + 1):
            matrizMinCuad[i][j] = np.dot(matrizPhis[i], matrizPhis[j])
            matrizMinCuad[j][i] = matrizMinCuad[i][j]

        terIndMinCuad[i] = np.dot(xtRuidoArr - gammaArr, matrizPhis[i])

    # Resolver para hallar valores de los coeficientes:
    return np.linalg.solve(matrizMinCuad, terIndMinCuad)    

# Reemplaza las expresiones simbolicas de los coeficientes por escalares y retorna lambdas de x(t) para todos los intervalos
def lambdificar(coefsResueltos):
    
    for cd in coefsDespejados:
        j = 0
        for ci in coefsIncognitas:
            ecsContYDeriv[cd] = ecsContYDeriv[cd].subs(ci, coefsResueltos[j])
            j += 1

    j = 0
    for ci in coefsIncognitas:
        ecsContYDeriv[ci] = coefsResueltos[j]
        j += 1

    global xLambdaArr
    for i in range(cantInts):
        xLambdaArr.append(sp.lambdify(t, xi(i)))

    return lambda t: xLambdaArr[obtenerIndiceIntervalo(t)](t)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Inciso 2

# Supone un pico con dos valles a los laterales. Las cotas deben contener a los dos minimos para asegurar el resultado correcto
def obtenerDuracionYAmplitud(piso, techo, tArr, xtAproxArr):
    
    # Quedarse solo con la porcion de la senial correspondiente a la onda:
    indPiso = indiceDelMasCercano(piso, tArr)
    indTecho = indiceDelMasCercano(techo, tArr)

    tArrOnda = tArr[indPiso : indTecho]
    xtAproxArrOnda = xtAproxArr[indPiso : indTecho] 
   
    # Minimos y maximos:
    indMinLocales = argrelmin(xtAproxArrOnda)[0]
    indMaxLocales = argrelmax(xtAproxArrOnda)[0]

    # Duracion:
    tMins = np.zeros(len(indMinLocales))
    for i in range(len(indMinLocales)):
        tMins[i] = tArrOnda[indMinLocales[i]]
    
    # Amplitud:
    xtMins = np.zeros(len(indMinLocales))
    for i in range(len(indMinLocales)):
        xtMins[i] = xtAproxArrOnda[indMinLocales[i]]
    
    xtMaxs = np.zeros(len(indMaxLocales))
    for i in range(len(indMaxLocales)):
        xtMaxs[i] = xtAproxArrOnda[indMaxLocales[i]]
    
    return tMins[len(tMins) - 1] - tMins[0], np.amax(xtMaxs) - np.amin(xtMins)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Inciso 3

# Calcula el ajuste minimo, promedio y maximo para cada intervalo, ademas del minimo, promedio y maximo general
def calcularAjuste(tArr, xtRuidoArr, xtAproxArr):

    # Ajustes para cada punto:
    ajustesArr = np.zeros(len(tArr))
    for i in range(len(tArr)):
        ajustesArr[i] = abs(xtAproxArr[i] - xtRuidoArr[i]) / (abs(xtAproxArr[i]) + 1) * 100   # https://stats.stackexchange.com/a/86710 

    # Indices divide-intervalos:
    # No tomar ni el primer ni el ultimo intervalo
    divInts = np.zeros(cantInts - 1, dtype = np.int64)
    for i in range (1, cantInts):
        divInts[i - 1] = np.where(extInts[i] == tArr)[0][0] # Retorna el indice del extremo

    # Para cada intervalo, calcular el minimo, promedio y maximo de los ajustes:
    ajustesInts = []
    for arr in np.split(ajustesArr, divInts):
        minAj = np.amin(arr)
        promAj = np.average(arr)
        maxAj = np.amax(arr)
        ajustesInts.append((minAj, promAj, maxAj))

    return ajustesInts, (np.amin(ajustesArr), np.average(ajustesArr), np.amax(ajustesArr))

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# General

def graficar(tArr, xtRuidoArr, xtAproxArr, xtLabel, title):
    plt.plot(tArr, xtRuidoArr, color = 'b', label = 'Con ruido')
    plt.plot(tArr, xtAproxArr, color = 'r', label = 'Mínimos cuadrados')
    plt.xlabel('t')
    plt.ylabel(xtLabel)
    plt.title(title)
    plt.legend()
    plt.show()

def aproximarSenial():
    
    #Inciso 1
    resolverEcsContYDeriv()    
    coefsResueltos = minimosCuadrados(tArr, xtRuidoArr)
    xLambda = lambdificar(coefsResueltos)
    xtAproxArr = np.array(list(map(xLambda, tArr)))

    graficar(tArr, xtRuidoArr, xtAproxArr, 'x(t)', 'ECG')
    
    print("x_minimosCuadrados(t) = ")
    for i in range(cantInts):
        print("\t{ " + str(xi(i)) + ",\t\t " + str(extInts[i]) + " <= t < " + str(extInts[i + 1]))

    #Inciso 2
    print("\n\nOnda\tDuración\t\tAmplitud")
    
    pisoP = 5. 
    techoP = 20.
    dur, amp = obtenerDuracionYAmplitud(pisoP, techoP, tArr, xtAproxArr)
    print("P\t" + str(dur) + "\t" + str(amp)) 

    pisoQRS = 15. 
    techoQRS = 30.
    dur, amp = obtenerDuracionYAmplitud(pisoQRS, techoQRS, tArr, xtAproxArr)
    print("QRS\t" + str(dur) + "\t" + str(amp))   

    pisoT = 22.
    techoT = 55.
    dur, amp = obtenerDuracionYAmplitud(pisoT, techoT, tArr, xtAproxArr)
    print("T\t" + str(dur) + "\t\t\t" + str(amp))  

    #Inciso 3
    ajustesInts, ajusteGen = calcularAjuste(tArr, xtRuidoArr, xtAproxArr)
    
    print("\n\n\t\t\t\tAjustes\n")
    print("Intervalo\tMínimo\t\t\tPromedio\t\tMáximo")
    i = 0
    for aj in ajustesInts:
        minAj, promAj, maxAj = aj
        print("[" + str(extInts[i]) + "; " + str(extInts[i + 1]) + ")\t" + str(minAj) + "\t" + str(promAj) + "\t" + str(maxAj))
        i += 1
    minAj, promAj, maxAj = ajusteGen
    print("\n[" + str(extInts[0]) + "; " + str(extInts[cantInts]) + ")\t" + str(minAj) + "\t" + str(promAj) + "\t" + str(maxAj) + "\n\n\n")

if __name__ == "__main__":
    aproximarSenial()