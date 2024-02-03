import sympy                 as sp
from   math import factorial as fact
from   os   import system    as sys

# ---------------------------------------------------------------------------------------------------------------------------------
x = sp.Symbol("x")
e = 1e-4
# xArr: valores de abscisa
# yArr: valores de ordenada; yArr[i] = p(xArr[i])
# dnyArr: [(x, [dy, d2y, d3y, ...])*]
def newtonInterpol(xArr, yArr, dnyArr):

    # Hallar derivadas a interpolar de una abscisa:
    def dnyDe(x):
        for j in range(0, len(dnyArr)):
            xj, dnyj = dnyArr[j]
            if abs(xj - x) < e:
                return dnyj
        return []

    # Repetir valores de los cuales se quieran interpolar derivadas:
    xArrMod = []
    yArrMod = []
    for j in range(0, len(xArr)):
        for k in range(0, 1 + len(dnyDe(xArr[j]))):
            xArrMod.append(xArr[j])
            yArrMod.append(yArr[j])

    # Armar en profundidad la tabla:
    difDiv = [yArrMod]
    i = 1
    # Mientras se puedan calcular diferencias divididas
    while(len(difDiv[i - 1]) >= 2):
        diy = []
        for j in range(0, len(xArrMod) - i):
            try:
                # Usar calculo normal, si no es 0/0
                m = (difDiv[i - 1][j + 1] - difDiv[i - 1][j]) / (xArrMod[j + i] - xArrMod[j])
            except:
                # Reemplazar con derivada / i!
                m = dnyDe(xArrMod[j])[i - 1] / fact(i)

            diy.append(m)
        
        # Anexar nivel i-esimo
        difDiv.append(diy)
        i += 1

    # Armar polinomio simbolico en reversa:
    p = 0
    for i in range(len(difDiv) - 1, -1, -1):
        p = difDiv[i][0] + (x - xArrMod[i]) * p
        
    return sp.expand(p)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
formatter = "{:4.4f}"
espaciosAdic = 5
def pEval(xArr, yArr, dnyArr, egNum):

    p = newtonInterpol(xArr, yArr, dnyArr) 
    print("p" + str(egNum) + "(x) = " + str(p) + "\n")

    for xj in xArr:
        print("p" + str(egNum) + "(" + str(xj) + ") = " + str(formatter.format(p.subs(x, xj))))

    print("")
    for xj, dny in dnyArr:
        for i in range(1, len(dny) + 1):
            print("p" + str(egNum) + ("\'" * i) + (" " * (espaciosAdic - i)) + "(" + str(formatter.format(xj)) + ") = " + str(formatter.format(p.diff(x, i).subs(x, xj))))
    
    # Supone que las abcisas vienen en orden creciente
    intAdic = 5.
    a = xArr[0] - intAdic
    b = xArr[len(xArr) - 1] + intAdic
    sp.plotting.plot(p, (x, a, b), xlabel = "x", ylabel = "p" + str(egNum) + "(x)")

    sys("cls")

if __name__ == "__main__":
    
    #Ejemplo 1
    xArr = [1., 2., 4., 6.]
    yArr = [1., 6., 6.5, 7.]
    dnyArr = []
    pEval(xArr, yArr, dnyArr, 1)

    #Ejemplo 2
    xArr = [1., 2., 4., 6.]
    yArr = [1., 4., 6.5, 7.]
    dnyArr = [(1., [9., 3., 7., 2.]), (2, [5., 7.])]
    pEval(xArr, yArr, dnyArr, 2)

    #Ejemplo 3
    xArr = [1.]
    yArr = [1.]
    dnyArr = [(1., [0.5, 1., 1.2])]
    pEval(xArr, yArr, dnyArr, 3)

    #Ejemplo 4
    xArr = [1., 2., 50.]
    yArr = [-1, 6., 6.5]
    dnyArr = [(1., [0.5, -0.8, 3.14, 8., 0.]), (2., [0.3, 1.])]
    pEval(xArr, yArr, dnyArr, 4)

    #Ejemplo 5
    xArr = [1., 2.]
    yArr = [1., 4.]
    dnyArr = [(1., [0.5, 1., 1.5])]
    pEval(xArr, yArr, dnyArr, 5)


