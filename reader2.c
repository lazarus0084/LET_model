#include <stdio.h>

int main() {
    FILE *fptr;
    char store[1000];

    fptr = fopen("input1.txt", "r");

    if (fptr != NULL) {
        while (fgets(store, 1000, fptr)) {
            printf("%s\n", store);
        }

        printf("file opened successfully");
    } else {
     printf("FILE NOT FOUND");
    }

    printf("\n\n\n");
    printf("completed\n");
    printf("checking...%s\n", store);

    return 0;
}
