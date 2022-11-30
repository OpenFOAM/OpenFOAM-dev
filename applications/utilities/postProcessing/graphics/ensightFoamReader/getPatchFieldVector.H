label patchi = which_part - 2;

if (nVar >= Num_variables - nSprayVariables)
{
    return Z_UNDEF;
}

typeIOobject<volVectorField> fieldObjectPtr
(
    fieldNames[var2field[nVar]],
    runTime.name(),
    mesh,
    IOobject::NO_READ
);

if (!fieldObjectPtr.headerOk())
{
    return Z_UNDEF;
}

IOobject fieldObject
(
    fieldNames[var2field[nVar]],
    runTime.name(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
);

volVectorField sf
(
    fieldObject,
    mesh
);

const vectorField& sfb = sf.boundaryField()[patchi];
const polyBoundaryMesh& bMesh = meshPtr->boundaryMesh();

if (which_type == Z_TRI03)
{
    label counter = 1;
    forAll(sfb, facei)
    {
        label nPoints = bMesh[patchi][facei].size();
        if (nPoints == 3)
        {
            var_array[counter++] = sfb[facei][component];
        }
    }
}

if (which_type == Z_QUA04)
{
    label counter = 1;
    forAll(sfb, facei)
    {
        label nPoints = bMesh[patchi][facei].size();
        if (nPoints == 4)
        {
            var_array[counter++] = sfb[facei][component];
        }
    }
}

if (which_type == Z_NSIDED)
{
    label counter = 1;
    forAll(sfb, facei)
    {
        label nPoints = bMesh[patchi][facei].size();
        if ((nPoints != 3) && (nPoints != 4))
        {
            var_array[counter++] = sfb[facei][component];
        }
    }
}
