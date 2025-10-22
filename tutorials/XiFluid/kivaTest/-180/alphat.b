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
    location    "-180";
    object      alphat.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [1 -1 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    piston
    {
        type            compressible::alphatWallFunction;
        Prt             0.85;
        value           uniform 0;
    }
    liner
    {
        type            compressible::alphatWallFunction;
        Prt             0.85;
        value           uniform 0;
    }
    cylinderHead
    {
        type            compressible::alphatWallFunction;
        Prt             0.85;
        value           uniform 0;
    }
}


// ************************************************************************* //
