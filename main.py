import imagini as img
import svd
import matplotlib.pyplot as plt
 
A = [[3, 1], [1, 3]]
k = 20
baze_cifre = {}
date_antrenare = img.incarca_date_mnist("mnist_train.csv",20)
for cifra in range(10):
    A_cifra = date_antrenare[cifra]
    U, sigma, V = svd.svd(A_cifra, k)

    # Reconstruieste prima imagine din A cu k componente
    reconstructie = [0.0] * 784
    for j in range(k):
        u_j = [U[i][j] for i in range(784)]
        v_0j = V[0][j]  # prima imagine, componenta j
        for i in range(784):
            reconstructie[i] += sigma[j] * u_j[i] * v_0j

    imagine = [reconstructie[r*28:(r+1)*28] for r in range(28)]
    plt.imshow(imagine, cmap="gray")
    plt.title(f"Cifra {cifra} reconstruita cu k={k}")
    plt.show()


print("S-a gatat codu sefule")