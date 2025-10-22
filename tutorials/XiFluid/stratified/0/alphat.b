/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  dev
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       volScalarField;
    location    "0";
    object      alphat.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [1 -1 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    left
    {
        type            symmetryPlane;
    }
    right
    {
        type            symmetryPlane;
    }
    top
    {
        type            symmetryPlane;
    }
    bottom
    {
        type            symmetryPlane;
    }
    frontAndBack
    {
        type            empty;
    }
}


// ************************************************************************* //
