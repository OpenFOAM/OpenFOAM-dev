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
    object      ft.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [];

internalField   uniform 0;

boundaryField
{
    #includeEtc "caseDicts/setConstraintTypes"

    inletFuel
    {
        type            inletOutlet;
        phi             alphaPhi.b;
        inletValue      uniform 0;
    }

    inlet
    {
        type            inletOutlet;
        phi             alphaPhi.b;
        inletValue      uniform 0;
    }

    outlet
    {
        type            inletOutlet;
        phi             alphaPhi.b;
        inletValue      uniform 0;
    }

    piston
    {
        type            zeroGradient;
    }

    liner
    {
        type            zeroGradient;
    }

    cylinderHead
    {
        type            zeroGradient;
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
