#include <iostream>
#include <unistd.h>

using namespace std;

int main()
{
    int *ptrs[500000];

//    for (;;);

    cerr << "allocating ints\n";

    for (int i=0; i<500000; i++)
    {
        ptrs[i] = new int[1];
        delete[] ptrs[i];
    }

    for (;;);

    cerr << "allocating double\n";

    double* array = new double[500000];

    for (;;);
}
