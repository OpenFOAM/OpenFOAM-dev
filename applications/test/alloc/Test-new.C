#include <stream.h>

main()
{
    int* intPtrs[500000];

    cerr << "allocating ints\n";

    for (int i=0; i<500000; i++)
    {
        intPtrs[i] = new int[1];
    }

    cerr << "allocated ints\n";

    cerr << "deallocating ints\n";

    for (i=0; i<500000; i++)
    {
        delete[] intPtrs[i];
    }

    cerr << "deallocated ints\n";

    cerr << "allocating doubles\n";

    double* doubles = new double[500000];

    for (;;);
}
