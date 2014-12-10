#include "stream.h"
#include <unistd.h>
#include <stdlib.h>

main()
{
    int *ptrs[500000];

    cerr << "allocating ints\n";

    for (int i=0; i<500000; i++)
    {
        ptrs[i] = (int*)malloc(sizeof(int));
    }

//    for (;;);

    cerr << "deallocating ints\n";

    for (i=0; i<500000; i++)
    {
        free(ptrs[i]);
    }

    cerr << "allocating double\n";

    double* array = (double*)malloc(500000*sizeof(double));

    for (;;);
}
