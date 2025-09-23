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
    class       dictionary;
    location    "constant";
    object      physicalProperties.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

thermoType
{
    type            ubPsiThermo;
    mixture         bInhomogeneousMixture;
    transport       polynomial;
    thermo          janaf;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleEnthalpy;
}

stoichiometricAirFuelMassRatio 17.23;

#include "$FOAM_CASE/constant/fuel_thermo.foam"

#include "$FOAM_CASE/constant/oxidant_thermo.foam"

#include "$FOAM_CASE/constant/products_thermo.foam"


// ************************************************************************* //
