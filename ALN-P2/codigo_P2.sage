def forma_escalonada(M, canon = False, ver_transfo = False, algoritmo = 'sin pivotear'):
    A  = copy(M)
    m  = A.nrows() # número de filas
    n  = A.ncols() # número de columnas
    mindim = min(m, n)
    lista_transfos = []                         # Lista que recordará las operaciones realizadas.
    F  = copy(identity_matrix(base_ring(A), m)) 
    s  = 0
    for k in range(mindim):
        (r, s) = busca_pivote[algoritmo](A,  k, s)
        if (r < m) and (s < n):   # Si hemos encontrado un pivote...
            if (r != k):          #    Si el pivote no está en la primera fila considerada...  
                A.swap_rows(k, r) #    ... lo colocamos allí (fila k)
                F.swap_rows(k, r)
                lista_transfos.append( ('swap',(k+1,r+1)) )
            piv = A[k, s]             # Valor del pivote. Su posicion es ahora: fila k, col s.
            ### reduccion debajo del pivote ### 
            for t in range(k + 1, m): # reducimos los coeficientes DEBAJO del pivote.
                multip = A[t, s]*piv^-1 
                if multip != 0:
                    A[t, s] = 0
                    F[t] = F[t] - multip*F[k]
                    lista_transfos.append( ('combination', (t+1, k+1, -multip)))
                    for v in range(s + 1, n):
                        A[t, v] = A[t, v] - multip*A[k, v]
            ### reduccion encima del pivote ###            
            if canon:    
                for t in range(k):   
                    multip = A[t, s]*piv^(-1)
                    if multip != 0:
                        A[t, s] = 0
                        F[t] = F[t] - multip*F[k]
                        lista_transfos.append( ('combination',(t+1, k+1, -multip)))
                        for v in range(s + 1, n):
                            A[t, v] = A[t, v] - multip*A[k, v]
                F[k] = F[k]*piv^(-1)  # reducimos también el pivote a 1.
                A[k] = A[k]*piv^(-1)
                lista_transfos.append( ('rescale', (k+1, piv^(-1))) )
        s = s+1        
    print
    ## Presentación de los resultados
    if ver_transfo: 
        show_operations(M, lista_transfos)
        if canon: print('Forma escalonada canónica:')
        else: print('Forma escalonada simple:')
        show(A)
    return(F, A)    
        
busca_pivote={}

def busca_pivote_sin_pivotear(A, k, s):
    '''Se busca como pivote el primer coeficiente no-nulo, 
    en las filas k y siguientes, columna s y siguientes.'''
    r = k
    (m, n) = A.dimensions()
    while (r < m) and (s < n) and (A[r, s] == 0):
       r = r + 1                    # ... en primer lugar en la misma columna s.
       if r == m: (r, s) = (k, s+1) # pasamos a la columna siguiente
    return(r,s)
    
busca_pivote['sin pivotear'] = busca_pivote_sin_pivotear

def show_operations(M, lista_transfos):
    'Presentacion de las operaciones elementales.'
    print('Matriz original:')
    show(M)
    print
    print('Transformaciones elementales realizadas:')
    for (tipo, parametros) in lista_transfos:
        if tipo == 'swap':
            (indice1, indice2) = parametros
            pretty_print(html('$F_{%s,%s}$'%(latex(indice1), latex(indice2))))
            print
        elif tipo == 'combination':
            (indice1, indice2, factor) = parametros
            pretty_print(html('$F_{%s,%s}(%s)$'%(latex(indice1),latex(indice2), latex(factor))))  
            print
        elif tipo == 'rescale':
            (indice, factor) = parametros
            pretty_print(html('$F_{%s}(%s)$'%(latex(indice), latex(factor))))  
            print
    return None

def sustitucion_regresiva(U):
    m = U.ncols()
    n = U.nrows()
    # Se separa la ultima columna: U = ( M | x )
    M = U.submatrix(ncols=m-1) # Las m-1 primeras columnas (se quita la ultima)
    x = U.column(-1)           # La ultima columna.
    for k in range(n - 1, -1, -1):
        x[k] = (x[k] - sum(M[k, j]*x[j] for j in range(k + 1, n)))/M[k, k]
    return x

def busca_pivote_con_pivoteo_parcial(A, k, s):
    vmax = 0    
    (m, n) = A.dimensions()
    r = m
    while (s < n) and (vmax == 0):
        # buscamos el coef max (en valor absoluto) en la col s, a partir de la fila k y por debajo
        (vmax, r) = max( (abs(A[p, s]), p)  for p in range(k, m) )  
        if vmax == 0: s = s + 1 # sin pivote en esta columna, pasamos a la siguiente.
    return (r, s)    
    
busca_pivote['pivoteo parcial'] = busca_pivote_con_pivoteo_parcial

def busca_pivote_con_pivoteo_parcial_escalado(A, k, s):
    vmax = 0    
    (m, n) = A.dimensions()
    r = m
    while (s < n) and (vmax == 0):
        for p in range(k,m):
            rowmax = max( abs(A[p,j]) for j in range(s,n) )
            if rowmax != 0 and abs(A[p,s])/rowmax > vmax: 
                (vmax, r) = (abs(A[p,s]) / rowmax , p)
        if vmax == 0: s = s+1
    return (r, s)    
    
    
busca_pivote['pivoteo parcial escalado'] = busca_pivote_con_pivoteo_parcial_escalado

def tril(A, k = 0):
    L = copy(A)
    m = A.nrows()
    n = A.ncols()
    L = matrix(A.base_ring(),m, n, lambda i,j : 0 if j > i+k else A[i,j])
    return(L)

def triu(A, k = 0):
    return(tril(A.transpose(), k).transpose())

def descompLU(M):
    A = copy(M)
    m = A.nrows(); n = A.ncols()
    for k in range(m):
        if A[k, k] != 0:
            for t in range(k + 1, m):
                A[t, k] = A[t, k]/A[k, k]
                A[t, k+1:] = A[t, k+1:] - A[t, k]*A[k, k+1:]
        elif any( x!=0 for x in A[k+1:, k]):
            print('La matriz no admite descomposición LU sin intercambio de filas')
            return None, None
    L = identity_matrix(base_ring(A), m) + tril(A, -1)
    U = triu(A)
    return L, U

def descompPLU(M):
    'Descomposicion PLU de A. Se tiene PA=LU.'
    A = copy(M)
    m = A.nrows(); n = A.ncols()
    piv = list(range(n))
    for k in range(m):
        (vmax, q) = max(  (abs(A[s, k]), s) for s in range(k, m) )
        if q != k:
            piv[k], piv[q] = piv[q], piv[k]
            A.swap_rows(q, k)
        if A[k, k] != 0:
            for t in range(k + 1, m):
                A[t, k]    = A[t, k]/A[k, k]
                A[t, k+1:] = A[t, k+1:] - A[t, k]*A[k, k+1:]
    L = identity_matrix(base_ring(A), m) + tril(A, -1)
    U = triu(A)
    P = identity_matrix(base_ring(A), m).matrix_from_rows(piv)
    return P, L, U
        
        
def descompPLUppe(M):
    A = copy(M)
    m = A.nrows(); n = A.ncols()
    piv = list(range(n))
    for k in range(m):
        vmax = 0
        for s in range(k,m):
            mm = max( abs(A[s,j]) for j in range(k,n) )
            if mm!=0:
                v = abs(A[s, k])/mm
                if v > vmax: (vmax, q) = (v, s)
        if q != k:
            piv[k], piv[q] = piv[q], piv[k]
            A.swap_rows(q, k)
        if A[k, k] != 0:
            for t in range(k + 1, m):
                A[t, k] = A[t, k]/A[k, k]
                A[t, k+1:] = A[t, k+1:] - A[t, k]*A[k, k+1:]
    L = identity_matrix(base_ring(A), m) + tril(A, -1)
    U = triu(A)
    P = identity_matrix(base_ring(A), m).matrix_from_rows(piv)
    return P, L, U

def modifica_lado_derecho(L, b):
    g = b.column()
    n = g.nrows()
    for k in range(n - 1):
        g[k+1:] = g[k+1:] - g[k,0]*L[k+1:, k]
    return g
    
        
def sustitucion_regresiva2(U, g):
    n = U.ncols()
    x = copy(g)
    for k in range(n - 1, -1, -1):
        x[k] = (g[k] - sum(U[k, j]*x[j] for j in range(k + 1, n)))/U[k, k]
    return x

def Cholesky(M):
    A = copy(M)
    R = A.base_ring()
    n = A.nrows()
    L = matrix(R, n, n, 0)
    if not A.is_symmetric():
        print('La matriz no es simétrica.')
        return None
    else:
        for i in range(n):
            for j in range(i):
                L[i,j]=1/L[j,j] * (A[i,j] - add(L[i,k]*L[j,k] for k in range(j)))
            aux = A[i,i] - sum(L[i,k]**2 for k in range(i))
            if aux <= 0:
                print('La matriz no es definida positiva. No admite descomposición de Cholesky.')
                return None
            L[i, i] = sqrt(aux)
        return L.transpose()

def sustitucion_progresiva2(L, g):
    n = L.ncols()
    x = copy(g)
    for k in range(n):
        x[k] = (g[k] - sum(L[k, j]*x[j] for j in range(k)))/L[k, k]
    return x



