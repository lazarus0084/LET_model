#include <stdio.h>

void save_my_mat(int arr[][3], int rows, const char *filename) {
    FILE *fp;
    fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Error opening file!\n");
        return;
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < 3; j++) {
            fprintf(fp, "%d ", arr[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    printf("Array saved to the file\n");
}

int main() {
    int A[][3] = {
        {1, 2, 3},
        {4, 5, 6},
        {8, 8, 8}
    };
    int size = sizeof(A) / sizeof(A[0]);
    save_my_mat(A, size, "gen_mat.txt");
    return 0;
}
