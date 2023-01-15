    // float *final_averages = malloc(file_count * NLAT * NLON * sizeof(float));
#include <stdio.h>
#include <stdlib.h>

int main()
{

	// This pointer will hold the
	// base address of the block created
	float* ptr;
	int n;

	// Get the number of elements for the array
	printf("Enter number of elements:");
	scanf("%d",&n);
	printf("Entered number of elements: %d\n", n);

	// Dynamically allocate memory using malloc()
	// ptr = (float*)malloc(n * sizeof(float));
    ptr = (float*)calloc(n, sizeof(float));

	// Check if the memory has been successfully
	// allocated by malloc or not
	if (ptr == NULL) {
		printf("Memory not allocated.\n");
		exit(0);
	}
	else {

		// Memory has been successfully allocated
		printf("Memory successfully allocated using malloc.\n");
        printf("%f\n", ptr[1]);
        // Get the elements of the array
        // for (i = 0; i < n; ++i) {
        // 	ptr[i] = i + 1;
        // }

        // Print the elements of the array
		printf("The elements of the array are: ");
		// for (i = 0; i < n; ++i) {
		// 	printf("%d, ", ptr[i]);
		// }
	}

	return 0;
}
