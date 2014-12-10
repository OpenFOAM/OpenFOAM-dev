scalar DiNum = -GREAT;

forAll(solidRegions, i)
{
    //- Note: do not use setRegionSolidFields.H to avoid double registering Cp
    //#include "setRegionSolidFields.H"
    const solidThermo& thermo = thermos[i];

    tmp<volScalarField> magKappa;
    if (thermo.isotropic())
    {
        magKappa = thermo.kappa();
    }
    else
    {
        magKappa = mag(thermo.Kappa());
    }

    tmp<volScalarField> tcp = thermo.Cp();
    const volScalarField& cp = tcp();

    tmp<volScalarField> trho = thermo.rho();
    const volScalarField& rho = trho();

    DiNum = max
    (
        solidRegionDiffNo
        (
            solidRegions[i],
            runTime,
            rho*cp,
            magKappa()
        ),
        DiNum
    );

}
