
def Fij(M, i, j):
    M.swap_rows(i - 1, j - 1)
    return M
    
def Fia(M, i, a):
    M.rescale_row(i - 1, a)
    return M
    
def Fija(M, i, j, a):
    M.add_multiple_of_row(i - 1, j -1, a)
    return M
    
def Cij(M, i, j):
    (M.transpose()).swap_rows(i - 1, j - 1)
    return M
    
def Cia(M, i, a):
    M.rescale_col(i - 1, a)
    return M
    
def Cija(M, i, j, a):
    M.add_multiple_of_column(i - 1, j - 1, a)
    return M

def forma_escalonada(M, canon = False, ver_transfo = False):
    """forma escalonada de una matriz - Práctica 1 de ALN."""
    A  = copy(M)
    m  = A.nrows() # número de filas
    n  = A.ncols() # número de columnas
    mindim = min(m, n)
    lista_transfos = []                         # Lista que recordará las operaciones realizadas.
    F  = copy(identity_matrix(base_ring(A), m)) 
    s  = -1
    for k in range(mindim):
        r = k      # r: indice de fila
        s = s + 1  # s: indice de columna
        # Buscamos pivote  en filas a partir de k ...
        while (r < m) and (s < n) and (A[r, s] == 0):
            r = r + 1             # ... en primer lugar en la misma columna s.
            if r == m:            # r ha superado el indice máximo de fila, que es m-1
                s = s + 1;r = k   # pasamos a la columna siguiente
        if (r < m) and (s < n):   # Si hemos encontrado un pivote...
            if (r != k):          #    Si el pivote no está en la primera fila considerada...  
                A.swap_rows(k, r) #    ... lo colocamos allí (fila k)
                F.swap_rows(k, r)
                lista_transfos.append( ('swap',(k+1,r+1)) )
            piv = A[k, s]             # Valor del pivote. Su posicion es ahora: fila k, col s.
            for t in range(k + 1, m): # reducimos los coeficientes DEBAJO del pivote.
                multip = A[t, s]*piv^-1 
                if multip != 0:
                    A[t, s] = 0
                    F[t] = F[t] - multip*F[k]
                    lista_transfos.append( ('combination', (t+1, k+1, -multip)))
                    for v in range(s + 1, n):
                        A[t, v] = A[t, v] - multip*A[k, v]
            if canon:    ## si queremos la forma canonica...
                for t in range(k):   # reducimos los coeficientes que estan ENCIMA del pivote
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
    print()
    ## Presentación de los resultados
    if ver_transfo:
        print('Matriz original:')
        print(M)
        print()
        print('Transformaciones elementales realizadas:')
        for (tipo, parametros) in lista_transfos:
            if tipo == 'swap':
                (indice1, indice2) = parametros
                #pretty_print(html('$F_{%s,%s}$'%(latex(indice1), latex(indice2))))
                pretty_print(html('$F_{%s} \\leftrightarrow F_{%s}$'%(latex(indice1), latex(indice2))))
                print()
            elif tipo == 'combination':
                (indice1, indice2, factor) = parametros
                #pretty_print(html('$F_{%s,%s}(%s)$'%(latex(indice1),latex(indice2), latex(factor))))
                pretty_print(html('$F_{%s}\\leftarrow F_{%s} + (%s) \\cdot F_{%s}$'
                                  %(latex(indice1),latex(indice1), latex(factor), latex(indice2))))
                print()
            elif tipo == 'rescale':
                (indice, factor) = parametros
                #pretty_print(html('$F_{%s}(%s)$'%(latex(indice), latex(factor))))
                pretty_print(html('$F_{%s} \\leftarrow %s\\cdot F_{%s}$'
                                  %(latex(indice), latex(factor), latex(indice))))
                print()
    if canon:
        print('Forma escalonada canónica:')
    else:
        print('Forma escalonada simple:')
    print(A)
    return(F, A)

