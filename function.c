#include <stdio.h>

int main() {
    // Declare the filename
   char *filename = "emptyyy.txt";

    // Open the file in write mode
    FILE *file = fopen(filename, "w");


    // Close the file
    fclose(file);

    printf("Empty file '%s' created successfully.\n", filename);

    return 0; // Exit the program successfully
}
