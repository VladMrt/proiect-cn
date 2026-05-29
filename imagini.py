
def incarca_date_mnist(path_fisier, limita_imagini_per_cifra=50):
    # Dictionar pentru a stoca matricea A pentru fiecare cifră
    # Cheia este cifra, valoarea este o lista de coloane
    matrici_clase = {i: [] for i in range(10)}
    contoare = {i: 0 for i in range(10)}

    with open(path_fisier, 'r') as f:
        # Trecem peste antet daca fisierul CSV are unul
        next(f, None) 
        
        for linie in f:
            valori = linie.strip().split(',')
            if not valori or len(valori) < 785:
                print("not valori sau len(valori)<785\n")
                continue
                
            cifra = int(valori[0]) # Primul element este eticheta (0-9)
            
            if contoare[cifra] >= limita_imagini_per_cifra:
                continue
                
            # Extragem pixelii si ii normalizam
            vector_imagine = [float(p) / 255.0 for p in valori[1:]]
            
            matrici_clase[cifra].append(vector_imagine)
            contoare[cifra] += 1
            
            if all(c >= limita_imagini_per_cifra for c in contoare.values()):
                break

    #facem dintr-o matrice 50 x 784 intr-una de 784 x 50
    for cifra in matrici_clase:
        matrice_pe_coloane = matrici_clase[cifra]
        n_coloane = len(matrice_pe_coloane)
        n_randuri = len(matrice_pe_coloane[0])
        
        A_final = [[matrice_pe_coloane[j][i] for j in range(n_coloane)] for i in range(n_randuri)]
        matrici_clase[cifra] = A_final

    return matrici_clase