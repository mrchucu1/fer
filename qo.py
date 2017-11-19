from sympy import Matrix, sympify
from utils.strings import sE
from utils.matrix import isQOrthogonal, illegalMatrix

#definiendo las funciones

def direccionesConjugadas(fStr, vars, Q, b, x0, ds):
    #le damos el formato que necesitamos operar
    Q = Matrix(Q)
    #validamos que nuestra Q y las direcciones que se nos pasan son ortogonales
    if not isQOrthogonal(ds, Q):
        raise illegalMatrix("La matriz no es simétrica.")
    b         = Matrix(b)
    x0        = Matrix(x0)
    iteracion = 0
    print("Usando el método de direcciones conjugadas", "con Q:\n%s\nb=%s\nd=\n%s\n en la siguiente función\n%s"%(sE(Q), sE(b.T), sE(x0.t), sE(ds), sE(fStr)))
    f = sympify(fStr)
    #misntras aún nos queden direcciones que operar iniciamos la iteración
    while len(ds) > 0:
        #obtenemos la dn
        d0 = ds.pop(0)
        iteracion += 1
        #generamos gn
        g0 = Q * x0 - b
        #generamos an
        a0 = -(g0.T*d0)/(d0.T*Q*d0)

        #obtenemos xn+1 apartir de xn y an
        x1    = x0 + a0[0] * d0
        zipX1 = list(zip(vars, x1))
        print("iteración %d:{\nx%d:\n%s\nf(x%d)=%f\n}\n"%(iteracion, iteracion, sE(x1), iteracion, f.subs(zipX1)))
        x0 = x1
    print("No. de iteraciones: ",iteracion)
    return zipX1
