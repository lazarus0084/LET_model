#include <stdio.h>

int main() {
    FILE* filePointer; // Declare a bookmark for our file
    filePointer = fopen("poem.txt", "r"); // Open the file called "poem.txt" for reading

    if (filePointer == NULL) {
        printf("Error: Can't open the file!\n");
        return 1; // If the file can't be opened, print an error message
    }

    // Now you can read from the file using filePointer
    char c;
    while ((c = fgetc(filePointer)) != EOF) { // Read each character until end of file
        putchar(c); // Print each character to the console
    }

    fclose(filePointer); // Close the file when you’re done
    return 0;
}
