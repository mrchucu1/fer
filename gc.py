from sympy import Matrix, sympify
from utils.strings import sE

#Definamos esta función para calcular el método del gradiente conjugado

def gradienteConjugado(fStr, vars, Q, b, x0):
    #Ponemos los datos a un formato que Python reconozca
    Q         = Matrix(Q)
    b         = Matrix(b)
    x0        = Matrix(x0)
    f0        = sympify(fStr)
    g0        = Q * x0 - b
    d0        = -1 * g0
    iteracion = 0
    print("Usando el método de gradiente conjugado ", "con Q:\n%s\nb=%s\nd0=\n%s\n en la siguiente función\n%s"%(sE(Q), sE(b.T), sE(x0.T), sE(d0.T), sE(fStr)))
    #Iteramos el número de variables que tengamos
    while true:
        iteracion += 1
        #Calculamos la an
        a0 = -(g0.T * d0)/(d0.T * Q * d0)
        #Calculamos xn+1 apartir de xn y an
        x1 = x0 + a0[0]*d0
        #Esta función es para asignar valores para el vector x1
        zipX1 = list(zip(vars, x1))
        print("Iteración %d:{\nx%d:\n%s\nf(x%d)=%s\n}\n"%(iteracion, iteracion, sE(x1), iteracion, f0.subs(zipX1)))
        #Verificamos si ya hicimos las iteraciones
        if iteracion == len(vars):
            #si ya salimos del ciclo
            break
        #De lo contrario si no hemos terminado con las iteraciones Calculamos
        # lo siguiente
        #calculamos gn+1
        g0 = Q * x1 - b
        #esto es auxiliar
        x0 = x1
        #calculamos la beta con la dirección dn
        beta = (g0.T * Q * d0) / (d0.T * Q * d0)
        #calculamos dn+1
        d0 = -1 * g0 + beta[0] * d0
    print("No. de iteraciones: ",iteracion)
    #al final regresamos el ultimo vector x calculado
    return zipX1
