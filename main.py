def transpusa(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]


def inmultire_matrici(A, B):
    """ A (m x n) @ B (n x p)"""
    m, n, p = len(A), len(A[0]), len(B[0])
    rezultat = [[0.0] * p for _ in range(m)] #rezultatul e o matrice m x p
    for i in range(m):
        for j in range(p):
            for k in range(n):
                rezultat[i][j] += A[i][k] * B[k][j]
    return rezultat

def modul(x):
    if x<0:
        return -x
    return x

def radical(x):
    return x ** 0.5


def matrice_identitate(n):
    I = [[0.0] * n for _ in range(n)]
    for i in range(n):
        I[i][i] = 1.0

def produs_scalar(v1, v2):
    """(dot product)"""
    suma = 0
    for i in range(len(v1)):
        suma += v1[i]*v2[i]
    return suma

#nu putem folosi good ol' iacobi de la algebra si geom sem1 ca matricea trebe sa fie simetrica :(
def descompunere_qr(A):
    """A = Q (ortogonala) x R (superior triunghiulara)."""
    n = len(A)

    Q = [[0.0] * n for _ in range(n)]
    R = [[0.0] * n for _ in range(n)]

   #modificam Q astfel incat sa devina ortogonala
    At = transpusa(A)
    Qt = [[0.0] * n for _ in range(n)]

    for i in range(n):
        v = list(At[i])#copie la coloana curenta

        for j in range(i):
            R[j][i] = produs_scalar(Qt[j], At[i])
            for k in range(n):
                v[k] -= R[j][i] * Qt[j][k]

        #R[i][i] = norma_vector(v) am sa o fac

        #evitam impartirea la 0
        if R[i][i] > 1e-12:
            for k in range(n):
                Qt[i][k] = v[k] / R[i][i]
        else:
            Qt[i][k] = 0.0

    Q = transpusa(Qt)
    return Q, R

A = [[1, 2], [3, 4], [5, 6]]
B = [[7],[8]]
print(transpusa(A))
print(modul(5))
print(modul(-7))
print(radical(5))
print(inmultire_matrici(A,B))
