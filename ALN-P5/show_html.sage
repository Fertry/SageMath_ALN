
def show_html(s):
    r"""
    Sustituye 'pretty_print(html(s))', 
    corrigiendo problemas que surgen con signos diacríticos 
    ÁÉÍÓÚÑáéíóúñ¿¡º
    en versiones recientes de Sage (por lo menos 8.2 y siguientes)
    """
    acentos = [('Á','A'),('É','E'),('Í','I'),('Ó','O'),('Ú','U'),('Ñ','N'),
               ('á','a'),('é','e'),('í','i'),('ó','o'),('ú','u'),('ñ','n')]
    for (XX, X) in acentos:
        s = s.replace(XX, '&{}acute;'.format(X)  )
    s = s.replace('¿','&iquest;').replace('¡','&iexcl;').replace('º','&ordm;').replace('−','-')
    pretty_print(html(s))

def matrix_from_copypaste(R, Nrows, Ncols, s):
    r"""
    Matriz del cuestionario.
    
    INPUT:
    - R -- anillo de base, por ejemplo RR, QQ, ZZ
    - Nrows -- numero de filas 
    - Ncols -- numero de columnas
    - s -- cadena de caracteres obtenida al copiar y pegar con el raton los coeficientes de un matriz del cuestionario
    """
    if len(s) > 0 and s[-1] == ",": 
        s = s[:-1]
    L = [ eval(x) for x in s.replace("−","-").split(",")]
    return matrix(R, Ncols, Nrows, L).transpose()
    

def vector_from_copypaste(R, s):
    r"""
    Vector del cuestionario.
    
    INPUT:
    - R -- anillo de base, por ejemplo RR, QQ, ZZ
    - s -- cadena de caracteres obtenida al copiar y pegar con el raton los coeficientes del vector del cuestionario
    """
    if len(s) > 0 and s[-1] == ",": 
        s = s[:-1]
    L = [ eval(x) for x in s.replace("−","-").split(",")]
    return vector(R, L)

