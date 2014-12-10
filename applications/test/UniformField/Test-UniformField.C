#include "UniformField.H"
#include "vector.H"
#include "IOstreams.H"

using namespace Foam;

int main()
{
    UniformField<scalar> uf1(13.1);
    UniformField<vector> uf2(vector(1, 2, 3));

    Info<< "uf1 = " << uf1[22] << "; uf2 = " << uf2[1002] << endl;

    return 0;
}
