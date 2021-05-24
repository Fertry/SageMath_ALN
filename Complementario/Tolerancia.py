################################################
# MODIFICAR SEGÚN LA REQUERIDA:
tolerancia = 2.9
nombre_fichero = 'imagen9.png'

################################################

import pylab
imagen = pylab.imread(DATA + nombre_fichero)
imagenmedia = pylab.mean(imagen, 2)

A = matrix(imagenmedia)
m = A.nrows()
n = A.ncols()
print("Matriz A de orden: ", m, 'x', n)

k = 0
norma = tolerancia + 1
while norma > tolerancia:
    k = k + 1
    Ak = U[:,0:k] * D[0:k,0:k] * (V[:,0:k].transpose())
    norma = (A - Ak).norm(2)

print("Nivel de tolerancia: ", tolerancia)
#compresion_imagen_SVD(k)
