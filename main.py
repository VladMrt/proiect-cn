import imagini as img
import svd
import matplotlib.pyplot as plt
import time

def proiecteaza(imagine_vector, U_baza):
    """
    Proiecteaza un vector imagine pe spatiul coloanelor lui U_baza.
    Returneaza distanta dintre imaginea originala si proiectia ei.
    imagine_vector: lista de 784 elemente
    U_baza: matrice 784 x k
    """
    m = len(U_baza)      # 784
    k = len(U_baza[0])   # numarul de vectori din baza

    # Calculam coeficientii proiectiei: c[j] = U[:,j] · imagine
    coef = [sum(U_baza[i][j] * imagine_vector[i] for i in range(m)) for j in range(k)]

    # Proiectia: p = U * c
    proiectie = [sum(U_baza[i][j] * coef[j] for j in range(k)) for i in range(m)]

    # Distanta = norma(imagine - proiectie)
    distanta = sum((imagine_vector[i] - proiectie[i]) ** 2 for i in range(m)) ** 0.5
    return distanta


# ── 1. Antrenare ──────────────────────────────────────────────────────────────
print("Antrenare SVD...")
start = time.time()

k = 10
baze_cifre = {}
date_antrenare = img.incarca_date_mnist("mnist_train.csv", 30)

for cifra in range(10):
    A_cifra = date_antrenare[cifra]
    U, sigma, V = svd.svd(A_cifra, k)
    baze_cifre[cifra] = U   # 784 x k

print(f"Antrenare terminata in {time.time() - start:.1f}s")


# ── 2. Testare ────────────────────────────────────────────────────────────────
primele_n = img.citeste_random_n("mnist_test.csv", n=25)
fig, axes = plt.subplots(5, 5, figsize=(8, 8))
corecte = 0

for idx, (cifra_reala, imagine_vec) in enumerate(primele_n):
    distante = {c: proiecteaza(imagine_vec, baze_cifre[c]) for c in range(10)}
    cifra_prezisa = min(distante, key=distante.get)
    culoare = "green" if cifra_prezisa == cifra_reala else "red"
    if cifra_prezisa == cifra_reala:
        corecte += 1

    # Convertim distantele in scoruri: scor = 1 - dist/sum(dist)
    suma_distante = sum(distante.values())
    scoruri = {c: (1 - distante[c] / suma_distante) * 100 for c in range(10)}

    imagine_2d = [imagine_vec[r*28:(r+1)*28] for r in range(28)]
    ax = axes[idx // 5][idx % 5]
    ax.imshow(imagine_2d, cmap="gray")
    ax.axis("off")
    ax.set_title(f"R:{cifra_reala} P:{cifra_prezisa}\n{scoruri[cifra_prezisa]:.1f}%",
                 fontsize=7, color=culoare)

plt.suptitle(f"Acuratete: {corecte}/25 ({100*corecte/25:.1f}%)", fontsize=12)
plt.tight_layout()
plt.show()