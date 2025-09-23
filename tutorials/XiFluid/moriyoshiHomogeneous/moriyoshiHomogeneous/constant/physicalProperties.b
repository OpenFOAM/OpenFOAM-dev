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
    object      physicalProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

thermoType
{
    type            ubPsiThermo;
    mixture         bHomogeneousMixture;
    transport       const;
    thermo          janaf;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleEnthalpy;
}

products
{
    specie
    {
        molWeight       28.3233;
    }
    thermodynamics
    {
        Tlow            200;
        Thigh           6000;
        Tcommon         1000;
        highCpCoeffs    (3.10561 0.00179748 -5.94701e-07 9.05612e-11 -5.08447e-15 -11003.6 5.12109);
        lowCpCoeffs     (3.498 0.000638554 -1.83885e-07 1.20991e-09 -7.68702e-13 -11080.6 3.1819);
    }
    transport
    {
        mu              1e-5;
        Pr              1;
    }
}

// ************************************************************************* //
