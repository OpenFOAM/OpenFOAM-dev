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
    object      T.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 1 0 0 0];

internalField   uniform 700;

boundaryField
{
    #includeEtc "caseDicts/setConstraintTypes"

    inletFuel
    {
        type            inletOutlet;
        phi             phic;
        inletValue      uniform 320;
        value           uniform 320;
    }

    inlet
    {
        type            inletOutlet;
        phi             phic;
        inletValue      uniform 350;
        value           uniform 350;
    }

    outlet
    {
        type            inletOutlet;
        phi             phic;
        inletValue      uniform 700;
        value           uniform 700;
    }

    piston
    {
        type            fixedValue;
        value           uniform 500;
    }

    liner
    {
        type            fixedValue;
        value           uniform 500;
    }

    cylinderHead
    {
        type            fixedValue;
        value           uniform 550;
    }

    ivHead
    {
        type            zeroGradient;
    }

    ivStem
    {
        type            zeroGradient;
    }

    evHead
    {
        type            zeroGradient;
    }

    evStem
    {
        type            zeroGradient;
    }

    frontAndBack
    {
        type            empty;
    }

    "nonCouple.*"
    {
        type            zeroGradient;
    }
}


// ************************************************************************* //
