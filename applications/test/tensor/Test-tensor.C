#include "tensor.H"
#include "symmTensor.H"
#include "transform.H"
#include "stringList.H"
#include "IOstreams.H"

using namespace Foam;

int main()
{
    tensor t1(1, 2, 3, 4, 5, 6, 7, 8, 9);
    tensor t2(1, 2, 3, 1, 2, 3, 1, 2, 3);

    tensor t3 = t1 + t2;

    Info<< t3 << endl;

    tensor t4(3,-2,1,-2,2,0,1, 0, 4);

    Info<< inv(t4) << endl;
    Info<< (inv(t4) & t4) << endl;

    Info<< t1.x() << t1.y() << t1.z() << endl;

    tensor t6(1,0,-4,0,5,4,-4,4,3);
    //tensor t6(1,2,0,2,5,0,0,0,0);
    Info<< "tensor " << t6 << endl;
    vector e = eigenValues(t6);
    Info<< "eigenvalues " << e << endl;
    tensor ev = eigenVectors(t6);
    Info<< "eigenvectors " << ev << endl;

    Info<< "Check determinant " << e.x()*e.y()*e.z() << " " << det(t6) << endl;

    Info<< "Check eigenvectors "
        << (eigenVector(t6, e[0]) & t6) << e[0]*eigenVector(t6, e[0]) << " "
        << (eigenVector(t6, e[1]) & t6) << e[1]*eigenVector(t6, e[1]) << " "
        << (eigenVector(t6, e[2]) & t6) << e[2]*eigenVector(t6, e[2])
        << endl;

    Info<< "Check eigenvalues for symmTensor "
        << eigenValues(symm(t6)) - eigenValues(tensor(symm(t6))) << endl;

    Info<< "Check eigenvectors for symmTensor "
        << eigenVectors(symm(t6)) - eigenVectors(tensor(symm(t6))) << endl;

    tensor t7(1, 2, 3, 2, 4, 5, 3, 5, 6);

    Info<< "Check transformation "
        << (t1 & t7 & t1.T()) << " " << transform(t1, t7) << endl;

    symmTensor st1(1, 2, 3, 4, 5, 6);
    symmTensor st2(7, 8, 9, 10, 11, 12);

    Info<< "Check symmetric transformation "
        << transform(t1, st1) << endl;

    Info<< "Check for dot product of symmetric tensors "
        << (st1 & st2) << endl;

    vector v1(1, 2, 3);

    Info<< sqr(v1) << endl;
    Info<< symm(t7) << endl;
    Info<< twoSymm(t7) << endl;
    Info<< magSqr(st1) << endl;
    Info<< mag(st1) << endl;

    Info<< (symm(t7) && t7) - (0.5*(t7 + t7.T()) && t7) << endl;
    Info<< (t7 && symm(t7)) - (t7 && 0.5*(t7 + t7.T())) << endl;


    /*
    // Lots of awkward eigenvector tests ...

    tensor T_rand_real
    (
        0.9999996423721313, 0.3330855667591095, 0.6646450161933899,
        0.9745196104049683, 0.0369445420801640, 0.0846728682518005,
        0.6474838852882385, 0.1617118716239929, 0.2041363865137100
    );
    Debug(T_rand_real);
    vector L_rand_real(eigenValues(T_rand_real));
    Debug(L_rand_real);
    tensor U_rand_real(eigenVectors(T_rand_real));
    Debug(U_rand_real);

    Info << endl << endl;

    tensor T_rand_imag
    (
        0.8668024539947510, 0.1664607226848602, 0.8925783634185791,
        0.9126510620117188, 0.7408077120780945, 0.1499115079641342,
        0.0936608463525772, 0.7615650296211243, 0.8953040242195129
    );
    Debug(T_rand_imag);
    vector L_rand_imag(eigenValues(T_rand_imag));
    Debug(L_rand_imag);
    tensor U_rand_imag(eigenVectors(T_rand_imag));
    Debug(U_rand_imag);

    Info << endl << endl;

    tensor T_rand_symm
    (
        1.9999992847442627, 1.3076051771640778, 1.3121289014816284,
        1.3076051771640778, 0.0738890841603279, 0.2463847398757935,
        1.3121289014816284, 0.2463847398757935, 0.4082727730274200
    );
    Debug(T_rand_symm);
    vector L_rand_symm(eigenValues(T_rand_symm));
    Debug(L_rand_symm);
    tensor U_rand_symm(eigenVectors(T_rand_symm));
    Debug(U_rand_symm);

    Info << endl << endl;

    symmTensor T_rand_Symm
    (
        1.9999992847442627, 1.3076051771640778, 1.3121289014816284,
                            0.0738890841603279, 0.2463847398757935,
                                                0.4082727730274200
    );
    Debug(T_rand_Symm);
    vector L_rand_Symm(eigenValues(T_rand_Symm));
    Debug(L_rand_Symm);
    tensor U_rand_Symm(eigenVectors(T_rand_Symm));
    Debug(U_rand_Symm);

    Info << endl << endl;

    tensor T_rand_diag
    (
        0.8668024539947510, 0, 0,
        0, 0.7408077120780945, 0,
        0, 0, 0.8953040242195129
    );
    Debug(T_rand_diag);
    vector L_rand_diag(eigenValues(T_rand_diag));
    Debug(L_rand_diag);
    tensor U_rand_diag(eigenVectors(T_rand_diag));
    Debug(U_rand_diag);

    Info << endl << endl;

    tensor T_repeated
    (
        0, 1, 1,
        1, 0, 1,
        1, 1, 0
    );
    Debug(T_repeated);
    vector L_repeated(eigenValues(T_repeated));
    Debug(L_repeated);
    tensor U_repeated(eigenVectors(T_repeated));
    Debug(U_repeated);

    Info << endl << endl;

    tensor T_repeated_zero
    (
        1, 1, 1,
        1, 1, 1,
        1, 1, 1
    );
    Debug(T_repeated_zero);
    vector L_repeated_zero(eigenValues(T_repeated_zero));
    Debug(L_repeated_zero);
    tensor U_repeated_zero(eigenVectors(T_repeated_zero));
    Debug(U_repeated_zero);

    Info << endl << endl;

    tensor T_triple
    (
        2, 0, 0,
        0, 2, 0,
        0, 0, 2
    );
    Debug(T_triple);
    vector L_triple(eigenValues(T_triple));
    Debug(L_triple);
    tensor U_triple(eigenVectors(T_triple));
    Debug(U_triple);
    */


    return 0;
}
