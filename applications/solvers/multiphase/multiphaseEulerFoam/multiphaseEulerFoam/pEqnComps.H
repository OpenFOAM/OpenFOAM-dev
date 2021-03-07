PtrList<fvScalarMatrix> pEqnComps(phases.size());

{
    PtrList<volScalarField> dmdts(fluid.dmdts());

    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];
        const volScalarField& alpha = phase;
        volScalarField& rho = phase.thermoRef().rho();

        pEqnComps.set(phasei, new fvScalarMatrix(p_rgh, dimVolume/dimTime));
        fvScalarMatrix& pEqnComp = pEqnComps[phasei];

        // Density variation
        if (!phase.isochoric() || !phase.pure())
        {
            pEqnComp +=
                (
                    fvc::ddt(alpha, rho) + fvc::div(phase.alphaRhoPhi())
                  - fvc::Sp(fvc::ddt(alpha) + fvc::div(phase.alphaPhi()), rho)
                )/rho;
        }

        // Mesh dilatation correction
        if (mesh.moving())
        {
            pEqnComp += fvc::div(mesh.phi())*alpha;
        }

        // Compressibility
        if (!phase.incompressible())
        {
            if (pimple.transonic())
            {
                const surfaceScalarField phid
                (
                    IOobject::groupName("phid", phase.name()),
                    fvc::interpolate(phase.thermo().psi())*phase.phi()
                );

                pEqnComp +=
                    correction
                    (
                        (alpha/rho)*
                        (
                            phase.thermo().psi()*fvm::ddt(p_rgh)
                          + fvm::div(phid, p_rgh)
                          - fvm::Sp(fvc::div(phid), p_rgh)
                        )
                    );

                pEqnComps[phasei].relax();
            }
            else
            {
                pEqnComp +=
                    (alpha*phase.thermo().psi()/rho)
                   *correction(fvm::ddt(p_rgh));
            }
        }

        // Option sources
        if (fvModels.addsSupToField(rho.name()))
        {
            pEqnComp -= (fvModels.source(alpha, rho) & rho)/rho;
        }

        // Mass transfer
        if (dmdts.set(phasei))
        {
            pEqnComp -= dmdts[phasei]/rho;
        }
    }
}
