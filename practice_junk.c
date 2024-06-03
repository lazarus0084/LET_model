#include <stdio.h>

int main() {
    int c, k, i, j;

    printf("Define number of columns: ");
    scanf("%d", &c);

    printf("Define number of rows: ");
    scanf("%d", &k);

    int A[k][c];

    for (i = 0; i < k; i++) {
        for (j = 0; j < c; j++) {
            printf("Enter element at position %d, %d: ", (i+1), (j+1));
            scanf("%d", &A[i][j]);
        }
    }

    printf("Entered matrix:\n");
    for (i = 0; i < k; i++) {
        for (j = 0; j < c; j++) {
            printf("%d ", A[i][j]);
        }
        printf("\n");
    }

    return 0;
}
