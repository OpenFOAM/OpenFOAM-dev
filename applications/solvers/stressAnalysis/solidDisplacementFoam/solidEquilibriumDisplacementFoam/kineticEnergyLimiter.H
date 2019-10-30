if (!(runTime.timeIndex() % 5))
{
    volScalarField kineticEnergy(magSqr(Dcorr));
    dimensionedScalar intKineticEnergy = fvc::domainIntegrate(kineticEnergy);
    Info<< "kineticEnergy = " << intKineticEnergy.value() << endl;

    volScalarField kineticPower((Dcorr - Dcorr.oldTime()) & Dcorr);
    dimensionedScalar intKineticPower = fvc::domainIntegrate(kineticPower);
    Info<< "kineticPower = " << intKineticPower.value() << endl;

    scalar smi = intKineticPower.value()/(intKineticEnergy.value() + vSmall);
    Info<< "smi = " << smi << endl;

    if (smi < -small)
    {
        Info<< "Resetting Dcorr to 0" << endl;
        Dcorr == dimensionedVector(Dcorr.dimensions(), Zero);
    }
}
