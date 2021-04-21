### Para la practica 3 de ALN.

##### forma_escalonada
## De la práctica 2

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