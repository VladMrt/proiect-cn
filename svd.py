def transpusa(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]

def modul(x):
    if x<0:
        return -x
    return x

def radical(x):
    return x ** 0.5

def inmultire_matrici(A, B):
    m, n, p = len(A), len(A[0]), len(B[0])
    rezultat = [[0.0] * p for _ in range(m)] #rezultatul e o matrice m x p
    for i in range(m):
        for j in range(p):
            for k in range(n):
                rezultat[i][j] += A[i][k] * B[k][j]
    return rezultat

def produs_scalar(v1, v2):
    suma = 0
    for i in range(len(v1)):
        suma += v1[i]*v2[i]
    return suma

def norma_vector(v):
    suma = 0.0
    for x in v:
        suma += x * x
    return radical(suma)

#nu putem folosi good ol' iacobi de la algebra si geom sem1 ca matricea trebe sa fie simetrica :(
def descompunere_qr(A):
    n = len(A)
    R = [[0.0] * n for _ in range(n)]
    At = transpusa(A)
    Qt = [[0.0] * n for _ in range(n)]
 
    for i in range(n):
        v = list(At[i])
 
        # scadem proiectiile pe vectorii ortonormali deja calculati
        for j in range(i):
            R[j][i] = produs_scalar(Qt[j], At[i])
            for k in range(n):
                v[k] -= R[j][i] * Qt[j][k]
 
        R[i][i] = norma_vector(v)
 
        # normalizam pentru a obtine vectorul ortonormal
        if R[i][i] > 1e-12:
            for k in range(n):
                Qt[i][k] = v[k] / R[i][i]
        else:
            for k in range(n):
                Qt[i][k] = 0.0
 
    Q = transpusa(Qt)
    return Q, R

# rezolva sistemul A*x = b prin eliminare Gaussiana cu pivot partial
def rezolva_sistem(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
 
    for col in range(n):
        # pivot partial
        max_row = max(range(col, n), key=lambda r: modul(M[r][col]))
        M[col], M[max_row] = M[max_row], M[col]
 
        pivot = M[col][col]
        if modul(pivot) < 1e-14:
            continue
 
        # Eliminam elementele de sub pivot
        for row in range(col + 1, n):
            factor = M[row][col] / pivot
            for j in range(col, n + 1):
                M[row][j] -= factor * M[col][j]
 
    # substitutie inversa (rezolvam de jos in sus)
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        if modul(M[i][i]) < 1e-14:
            x[i] = 0.0
        else:
            x[i] = M[i][n]
            for j in range(i + 1, n):
                x[i] -= M[i][j] * x[j]
            x[i] /= M[i][i]
    return x

# calculeaza valorile proprii ale matricei simetrice A prin iteratii QR
def valori_proprii_qr(A, iteratii=1000, tol=1e-10):
    Ak = [row[:] for row in A]
    n = len(A)
 
    for _ in range(iteratii):
        Q, R = descompunere_qr(Ak)
        Ak = inmultire_matrici(R, Q)

        suma_off = 0.0
        for i in range(1, n):
            for j in range(i):
                suma_off += Ak[i][j] ** 2
        if suma_off < tol:
            break

    return [Ak[i][i] for i in range(n)]

def calculeaza_vectori_proprii(A, valori_proprii, iteratii=1000, tol=1e-12):
    n = len(A)
    V = []
 
    for lam in valori_proprii:
        v = [1.0 / radical(n)] * n  # vector initial uniform normalizat
        shift = lam + 1e-9          # shift mic pentru a evita singularitatea exacta a lui (A - lam*I)
 
        for _ in range(iteratii):
            # Construim (A - shift*I) si rezolvam sistemul pentru a aplica inversa
            B = [[A[i][j] - (shift if i == j else 0.0) for j in range(n)] for i in range(n)]
            w = rezolva_sistem(B, v)
            norma = norma_vector(w)
            if norma < tol:
                break
            v_nou = [x / norma for x in w]
            # Criteriu de convergenta: vectorul nu se mai schimba semnificativ
            diff = sum((v_nou[k] - v[k]) ** 2 for k in range(n))
            v = v_nou
            if diff < tol:
                break
 
        # Gram-Schmidt
        for vj in V:
            coef = sum(v[k] * vj[k] for k in range(n))
            v = [v[k] - coef * vj[k] for k in range(n)]
        norma = norma_vector(v)
        v = [x / norma for x in v] if norma > tol else v
 
        V.append(v)
 
    # Transpunem: V[i][j] = componenta i a vectorului propriu j
    return [[V[j][i] for j in range(len(V))] for i in range(n)]

def svd(A):
    m, n = len(A), len(A[0])
    At = transpusa(A)
 
    # valorile singulare
    AtA = inmultire_matrici(At, A)
    val_proprii = valori_proprii_qr(AtA)
    val_proprii.sort(reverse=True)
    valori_singulare = [round(radical(modul(x)), 10) for x in val_proprii] 
    # coloanele lui V = vectorii proprii ai A^T * A
    V = calculeaza_vectori_proprii(AtA, val_proprii)
 
    # coloanele lui U = A*v_i / sigma_i
    U = []
    for i in range(len(valori_singulare)):
        if valori_singulare[i] > 1e-12:
            Av = [sum(A[r][c] * V[c][i] for c in range(n)) for r in range(m)]
            norma = norma_vector(Av)
            U.append([x / norma for x in Av])
        else:
            u = [0.0] * m
            u[len(U)] = 1.0
            U.append(u)
 
    U = transpusa(U)  # am construit U pe linii, il transpunem ca coloanele sa fie vectori
    return U, valori_singulare, V