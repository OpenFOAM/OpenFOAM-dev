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
    #includeEtc "caseDicts/setConstraintTypes"

    inletFuel
    {
        type            calculated;
        value           uniform 0;
        Prt             0.85;
    }

    inlet
    {
        type            calculated;
        value           uniform 0;
        Prt             0.85;
    }

    outlet
    {
        type            calculated;
        value           uniform 0;
        Prt             0.85;
    }

    piston
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
        Prt             0.85;
    }

    liner
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
        Prt             0.85;
    }

    cylinderHead
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
        Prt             0.85;
    }

    ivHead
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
        Prt             0.85;
    }

    ivStem
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
        Prt             0.85;
    }

    evHead
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
        Prt             0.85;
    }

    evStem
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
        Prt             0.85;
    }

    frontAndBack
    {
        type            empty;
    }

    "nonCouple.*"
    {
        type            calculated;
        value           uniform 0;
    }
}


// ************************************************************************* //
