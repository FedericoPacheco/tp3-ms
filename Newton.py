import math
import sympy as sp
from sympy import *
from math import factorial

X = Symbol('x')

# -------------------- Funciones -----------------

def tablaDiferenciasDivididas(m, x, derivadas):
  
    
    for i in range(1, len(m)): #columna
        for j in range(len(m) - i): #fila

            numerador = (m[j + 1][i - 1] - m[j][i - 1])
            divisor = (x[i + j] - x[j])

            if(divisor != 0):
                m[j][i] = ( numerador / divisor)
            else:
                m[j][i] = derivadas[i-1].get(x[j]) / factorial(i)

    
    return m

# -----------------------------------------------

def cantidadDeCondiciones(derivadas):
    n = 0
    for i in range(len(derivadas)):
        n += len(derivadas[i])
    return n

# -----------------------------------------------

def inicializarTabla(x, y, derivadas, n):
    
    tabla = [[0 for i in range(n+len(y))] 
                for j in range(n+len(y))]
    xActualizado = []
    k = 0
    i = 0
    while(i < n+len(y)):
    #for i in range (0,len(x)):
        for j in range (0, len(derivadas)):    
            if( x[k] in list(derivadas[j].keys()) ):
                tabla[i][0] = y[k]
                xActualizado.append(x[k])
                i+=1
        
        tabla[i][0] = y[k]
        xActualizado.append(x[k])
        k+=1
        i+=1

    return tabla, xActualizado

# -----------------------------------------------

def printTabla(tabla): 
  
    for i in range(len(tabla)): 
        for j in range(len(tabla) - i): 
            print(round(tabla[i][j], 4), "\t", 
                               end = " ")
  
        print("") 

# -----------------------------------------------

def armarPolinomio(x, tabla): 
    
    sum = tabla[0][0]; 
  
    for i in range(1, len(tabla)):
        sum = sum + (tabla[0][i] * proterm(i, x)); 

    return sum

# -----------------------------------------------
# Funcion para hallar el producto de los terminos 
def proterm(i, x): 
    global X
    pro = 1; 
    for j in range(i): 
        pro = pro * (X - x[j]); 
    
    return pro

# -----------------------------------------------

def interpolacionNewton(x, y, derivadas):

  # x es el vector de abscisas, y es el vector de ordenadas
  # derivadas es un vector de diccionarios, 
  # donde la clave es la abscisa y el valor es el valor de la derivada correspondiente a
  # la posicion del arreglo(+1). Se asume que si se busca interpolar la derivada segunda en 
  # un punto, tambien se quiere interpolar la primera.
    
    n = cantidadDeCondiciones(derivadas) #Condiciones adicionales

    m, xActualizado = inicializarTabla(x, y, derivadas, n)

    tabla = tablaDiferenciasDivididas(m, xActualizado, derivadas)
    print("Tabla de diferencias divididas:")
    printTabla(tabla)
    print("")

    print("Polinomio")
    ec = armarPolinomio(xActualizado, tabla)
    #print("P(x) = " + str(ec))
    print("P(x) = " + str(simplify(ec)))
    print("")
    
    print("Puntos:")
    for i in range(len(x)):
        print("P(" + str(x[i]) + ") = " + str(ec.subs(X,x[i])))
        if(i!=0 and x[i] == x[i-1]): i+=1
    
    print("")

    print("Derivadas:")
    derivada = ec.diff(X)
    for j in range(len(derivadas)):
        for k in range(len(derivadas[j])):
            abs = list(derivadas[j].keys())[k]
            ord = derivada.subs(X,abs)
            print("P^(" + str(j+1) + ")(" + str(abs) + ") = " + str(ord) )
        derivada = derivada.diff(X)

    print("")

    sp.plotting.plot(ec,(X,x[0]-1,x[len(x)-1]+1), title='Polinomio de Newton', xlabel='x', ylabel='P(x)')
    

# -----------------------------------------------

if __name__ == "__main__":
    
    x = [1, 2, 4, 6]
    y = [1, 6, 6.5, 7]

    # derivadas es un vector de diccionarios, 
    # donde la clave es la abscisa y el valor es el valor de la derivada correspondiente a
    # la posicion del arreglo(+1). 
    # Se asume que si se busca interpolar la derivada segunda en 
    # un punto, tambien se quiere interpolar la primera.

    #Ejemplo 1
    print("---------------------------Ejemplo 1---------------------------")
    deriv = [] #Sin restricciones de derivada
    interpolacionNewton(x,y, deriv)

    #Ejemplo 2
    print("---------------------------Ejemplo 2---------------------------")
    deriv = [ {1:9, 2:5}, # Esto corresponde a la primera derivada
              {1:3, 2:7}, # Esto corresponde a la segunda derivada
              {1:7},      # Esto corresponde a la tercera derivada
              {1:2}       # Esto corresponde a la cuarta derivada
            ]

    interpolacionNewton(x,y, deriv)

    #Ejemplo 3
    print("---------------------------Ejemplo 3---------------------------")
    x = [1]
    y = [1]
    deriv = [ {1:0.5}, {1:1}, {1:1.2} ]
    interpolacionNewton(x,y, deriv)

    #Ejemplo 4
    print("---------------------------Ejemplo 4---------------------------")
    x = [1, 2, 50]
    y = [-1, -6, -6.5]

    deriv = [ {1:0.5, 2:0.3}, # Esto corresponde a la primera derivada
              {1:-0.8, 2:1},  # Esto corresponde a la segunda derivada
              {1:math.pi},    # Esto corresponde a la tercera derivada
              {1:8},          # Esto corresponde a la cuarta derivada
              {1:0},          # Esto corresponde a la quinta derivada
            ]

    interpolacionNewton(x,y, deriv)


    #Ejemplo 5
    print("---------------------------Ejemplo 5---------------------------")
    x = [1, 2]
    y = [1, 4]
    deriv = [{1:0.5}, {1:1}, {1:1.5}]
    interpolacionNewton(x,y, deriv)