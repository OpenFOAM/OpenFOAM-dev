label patchi = which_part - 2;

if (nVar >= Num_variables - nSprayVariables)
{
    return Z_UNDEF;
}

IOobject fieldObjectPtr
(
    fieldNames[var2field[nVar]],
    runTime.timeName(),
    mesh,
    IOobject::NO_READ
);

if (!fieldObjectPtr.typeHeaderOk<volTensorField>(true))
{
    return Z_UNDEF;
}

IOobject fieldObject
(
    fieldNames[var2field[nVar]],
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
);

volTensorField sf
(
    fieldObject,
    mesh
);

const tensorField& tf = sf.boundaryField()[patchi];
const polyBoundaryMesh& bMesh = meshPtr->boundaryMesh();

if (which_type == Z_TRI03)
{
    label counter = 1;
    forAll(tf, n)
    {
        label nPoints = bMesh[patchi][n].size();
        if (nPoints == 3)
        {
            #include "tensorConversion.H"
        }
    }
}

if (which_type == Z_QUA04)
{
    label counter = 1;
    forAll(tf, n)
    {
        label nPoints = bMesh[patchi][n].size();
        if (nPoints == 4)
        {
            #include "tensorConversion.H"
        }
    }
}

if (which_type == Z_NSIDED)
{
    label counter = 1;
    forAll(tf, n)
    {
        label nPoints = bMesh[patchi][n].size();
        if ((nPoints != 3) && (nPoints != 4))
        {
            #include "tensorConversion.H"
        }
    }
}
