%auto
def matriz_relajacion(A, w):
    k = 1
    m = A.nrows()
    n = A.ncols()
    if m != n:
        k = 0
        print("La matriz no es cuadrada colega")
    else:
        for i in range(m):
            k = k*A[i, i]
        if k == 0:
            print("Es necesario permutar las filas de la matriz!")
        else:
            k = 1
            D = diagonal_matrix([A[j, j] for j in range(m))
            E = -tril(A, -1)
            F = -triu(A, -1)
            DwE1 = (D - w*E).inverse()
            return(DwE1*((1 - w)*D + w*F))


%auto
def matriz_gauss_seidel(A):
    k = 1
    m = A.nrows()
    n = A.ncols()
    if m != n:
        k = 0
        print("La matriz no es cuadrada colega")    
    else:
        for i in range(m):
            k = k*A[i, i]
        if k == 0:
            print("Es necesario permutar las filas de la matriz!!!")
        else:
            DE1 = tril(A, 0).inverse()
            F = -triu(A, -1)
            return(DE1*F)


%auto
def matriz_jacobi(A):
    k = 1
    m = A.nrows()
    n = A.ncols()
    if m !=  n:
        k = 0
        print("La matriz no es cuadrada colega")        
    else:
        for i in range(m):
            k = k*A[i, i]
        if k == 0:
            print("Es necesario permutar las filas de la matriz!!!")
        else:
            D = diagonal_matrix([A[j, j] for j in range(m)])
            D1 = diagonal_matrix([1/A[j, j] for j in range(m)])
            EF = -A + D
            J = D1*(EF)
            return(J)


%auto
def tril(A, k = 0):
    L = copy(A)
    m = L.nrows()
    n = L.ncols()
    for i in range(m):
        for j in range(i + 1 + k, n):
            L[i, j] = 0
    return(L)

def triu(A, k = 0):
    L = copy(A).transpose()
    return(tril(L, k).transpose())

def descomposicion(A):
    m = A.nrows()
    n = A.ncols()
    n = min(m, n)
    D = diagonal_matrix([A[i, i] for i in range(n)])
    E = -tril(A, -1)
    F = -triu(A, -1)


%auto
def radio_espectral(B):
    A = copy(B)
    return vector((A.change_ring(CC)).eigenvalues()).norm(Infinity)


%auto
def metodo_relajacion(A, b, w, tol = 10^-14, maxiter = 1000, veriter = False, solosiconverge = True):
    k = 1
    m = A.nrows()
    n = A.ncols()
    if m != n:
        k = 0
        print("La matriz no es cuadrada colega")
    else: 
        for i in range(m):
            k = k*A[i, i]
        if k == 0:
            print("Es necesario permutar las ecuaciones (no puede haber ceros en la diagonal)")
        else:
            k = 1; n = 0
            Lw = matriz_relajacion(A, w)
            R = radio_espectral(Lw).n()
            if solosiconverge and R >= 1:
                k = 0
                print("El método NO es convergente miarma: ")
                print("El radio espectral de la matriz de relajación es: ", R)
            else:
                x0 = vector(A.base_ring(), [0 for h in range(m)])
                x1 = copy(x0)
                while (A*x0 - b).norm() > tol and n < maxiter:
                    n = n + 1
                    for i in range(m):
                        si = 0; ss = 0
                        for t in range(i):
                            si = si + A[i, t]*x1[t]
                        for t in range(i + 1, m):
                            ss = ss + A[i, t]*x0[t]
                        x1[i] = w/A[i, i]*(b[i] - si - ss) + (1 - w)*x0[i]
                    x0 = copy(x1)
                    if veriter:
                        print("Iteración: ", n, " = ", x0)
    if k == 1:
        if n < maxiter:
            print("Nº de iteraciones realizadas: ", n) 
        else:
            print("Abortado en la iteración: ", maxiter)
        print("El radio espectral de la matriz de relajación es: ", R)
        Error = (A*x0 - b).norm().n()
        print("Error de la aproximación: ", Error)


%auto
def metodo_gauss_seidel(A, b, tol = 10^-14, maxiter = 1000, veriter = False, solosiconverge = True):
    k = 1
    m = A.nrows()
    n = A.ncols()
    if m != n:
        k = 0
        print("La matriz no es cuadrada colega")
    else: 
        for i in range(m):
            k = k*A[i, i]
        if k == 0:
            print("Es necesario permutar las ecuaciones (no puede haber ceros en la diagonal)")
        else:
            k = 1; n = 0
            L = matriz_gauss_seidel(A)
            R = radio_espectral(L).n()
            if solosiconverge and R >= 1:
                k = 0
                print("El método NO es convergente miarma: ")
                print("El radio espectral de la matriz de Gauss Seidel es: ", R)
            else:
                x0 = vector(A.base_ring(), [0 for h in range(m)])
                x1 = copy(x0)
                while (A*x0 - b).norm() > tol and n < maxiter:
                    n = n + 1
                    for i in range(m):
                        si = 0; ss = 0
                        for t in range(i):
                            si = si + A[i, t]*x1[t]
                        for t in range(i + 1, m):
                            ss = ss + A[i, t]*x0[t]
                        x1[i] = 1/A[i, j]*(b[i] - si - ss)
                    x0 = copy(x1)
                    if veriter:
                        print("Iteración: ", n, " = ", x0)
    if k == 1:
        if n < maxiter:
            print("Nº de iteraciones realizadas: ", n) 
        else:
            print("Abortado en la iteración: ", maxiter)
        print("El radio espectral de la matriz de Gauss Seidel es: ", R)
        Error = (A*x0 - b).norm().n()
        print("Error de la aproximación: ", Error)


%auto
def metodo_jacobi(A, b, tol = 10^-14, maxiter = 1000, veriter = False, solosiconverge = True):
    k = 1
    m = A.nrows()
    n = A.ncols()
    if m != n:
        k = 0
        print("La matriz no es cuadrada colega")
    else: 
        for i in range(m):
            k = k*A[i, i]
        if k == 0:
            print("Es necesario permutar las ecuaciones (no puede haber ceros en la diagonal)")
        else:
            k = 1; n = 0
            J = matriz_jacobi(A)
            R = radio_espectral(J).n()
            if solosiconverge and R >= 1:
                k = 0
                print("El método NO es convergente miarma: ")
                print("El radio espectral de la matriz de Jacobi es: ", R)
            else:
                x0 = vector(A.base_ring(), [0 for h in range(m)])
                x1 = copy(x0)
                while (A*x0 - b).norm() > tol and n < maxiter:
                    n = n + 1
                    for i in range(m):
                        si = 0; ss = 0
                        for t in range(i):
                            si = si + A[i, t]*x0[t]
                        for t in range(i + 1, m):
                            ss = ss + A[i, t]*x0[t]
                        x1[i] = 1/A[i, i]*(b[i] - si - ss)
                    x0 = copy(x1)
                    if veriter:
                        print("Iteración: ", n, " = ", x0)
    if k == 1:
        if n < maxiter:
            print("Nº de iteraciones realizadas: ", n) 
        else:
            print("Abortado en la iteración: ", maxiter)
        print("El radio espectral de la matriz de Jacobi es: ", R)
        Error = (A*x0 - b).norm().n()
        print("Error de la aproximación: ", Error)