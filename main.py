import imagini
import svd
 
A = [[3, 1], [1, 3]]
baze_cifre = {}
k = 10 
date_antrenare = imagini.incarca_date_mnist("/home/dutu/proiect-cn/mnist_test.csv",20)
for cifra in range(10):
    A_cifra = date_antrenare[cifra]
    U, sigma, V = svd.svd(A_cifra)
    # Pastram doar primele k coloane din U
    U_truncat = [[U[i][j] for j in range(k)] for i in range(len(U))]
    print(U_truncat)
    print("\n")
    baze_cifre[cifra] = U_truncat

print("S-a gatat codu sefule")