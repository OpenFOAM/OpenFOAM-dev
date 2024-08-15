
// Not sure if this is necessary anymore

nVar -= Num_variables - nSprayVariables + lagrangianScalarNames.size();

if (nVar >= 0)
{
    word name = lagrangianVectorNames[nVar];

    IOField<vector> v
    (
        IOobject
        (
            name,
            runTime.name(),
            lagrangian::cloud::prefix,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (v.size())
    {
        forAll(v, n)
        {
            if (component == 0)
            {
                var_array[n+1] = v[n].x();
            }
            else if (component == 1)
            {
                var_array[n+1] = v[n].y();
            }
            else if (component == 2)
            {
                var_array[n+1] = v[n].z();
            }
        }
    }
}
else
{
    // Info<< "getLagrangianVector: nVar = " << nVar << endl;
    return Z_UNDEF;
}
