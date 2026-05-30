import numpy as np
import matplotlib.pyplot as plt

date_mnist = np.loadtxt("mnist_train.csv", delimiter=",",skiprows=1)

etichete = date_mnist[:, 0].astype(int)#formatul MNIST contine la inceputul fiecarui rand o eticheta cu cifra careia ii apartin pixelii

date = date_mnist[:, 1:]

date = date.reshape(-1, 28, 28)
lungime_date = len(date)

imagini_originale = {}
imagini_reconstruite = {}
k = 10

for i in range(lungime_date):
    if len(imagini_originale) == 10:
        break
    cifra_actuala = etichete[i]
    if cifra_actuala not in imagini_originale:
        imagini_originale[cifra_actuala] = date[i]

        U, sigma, V_T = np.linalg.svd(imagini_originale[cifra_actuala], full_matrices=False)

        U_k = U[:, :k]
        sigma_k = np.diag(sigma[:k])
        V_T_k = V_T[:k, :]

        reconstructie = U_k @ sigma_k @ V_T_k
        imagini_reconstruite[cifra_actuala] = reconstructie

for cifra in range(10):

    plt.figure(figsize=(10, 4))
    plt.imshow(imagini_reconstruite[cifra], cmap="gray")
    plt.title(f"Reconstruita cu k={k}")

    plt.show()