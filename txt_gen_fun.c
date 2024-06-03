#include <stdio.h>

void save_my_mat(int arr[], int size, const char *filename) {
    FILE *fp;
    fp = fopen(filename, "w");

    for (int i = 0; i < size; i++) {
        fprintf(fp, "%d\n", arr[i]);
    }
    fclose(fp);
    printf("Array saved to the file\n");
}

int main() {
    int A[] = {1, 2, 3,4,6,5,6};
    int size = sizeof(A) / sizeof(A[0]);
    save_my_mat(A, size, "gen_arr.txt");
    return 0;
}
