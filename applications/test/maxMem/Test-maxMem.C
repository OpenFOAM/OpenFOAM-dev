#include <iostream>
#include <stdlib.h>
using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cerr << "Usage: " << argv[0] << " <number of Mb per chunk>\n";
        exit(1);
    }

    int nBytes = (1024U*1024U)*atoi(argv[1]);

    char *cPtr;

    for (unsigned i=1;; i++)
    {
        cPtr = new char[nBytes];

        /*
        for (int j=0; j<nBytes; j++)
        {
            cPtr[j] = 0;
        }
        */

        cout << "allocated " << i*nBytes/(1024U*1024U) << " Mbytes" << endl;
    }

    return 0;
}
