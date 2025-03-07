#include <stdio.h>
#include <stdlib.h>

void teste(int *m){
    for(int i = 0; i < 4; i++){
        m[i] = 0;
    }
}

int main() {
    int linhas = 3, colunas = 4;
    
    // Aloca um array de ponteiros para as linhas
    int **matriz = (int **)malloc(linhas * sizeof(int *));
    
    if (matriz == NULL) {
        printf("Erro ao alocar memória\n");
        return 1;
    }

    // Aloca um array para cada linha
    for (int i = 0; i < linhas; i++) {
        matriz[i] = (int *)malloc(colunas * sizeof(int));
        if (matriz[i] == NULL) {
            printf("Erro ao alocar memória para a linha %d\n", i);
            return 1;
        }
    }
    
    

    // Preenchendo a matriz
    for (int i = 0; i < linhas; i++) {
        for (int j = 0; j < colunas; j++) {
            matriz[i][j] = i * colunas + j;
        }
    }
    
    teste(matriz[0]);

    // Imprimindo a matriz
    for (int i = 0; i < linhas; i++) {
        for (int j = 0; j < colunas; j++) {
            printf("%d ", matriz[i][j]);
        }
        printf("\n");
    }

    // Liberando a memória
    for (int i = 0; i < linhas; i++) {
        free(matriz[i]);
    }
    free(matriz);
    
    return 0;
}
