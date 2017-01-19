
// Not sure if this is necessary anymore
nVar -= Num_variables - nSprayVariables;

if (nVar >= 0)
{
    word name = lagrangianScalarNames[nVar];

    IOField<scalar> s
    (
        IOobject
        (
            name,
            runTime.timeName(),
            cloud::prefix,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (s.size())
    {
        forAll(s, n)
        {
            var_array[n+1] = s[n];
        }
    }
}
else
{
    // Info<< "getLagrangianScalar: nVar = " << nVar << endl;
    return Z_UNDEF;
}
