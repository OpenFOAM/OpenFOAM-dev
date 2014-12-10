{
    // rho1 = rho10 + psi1*p;
    // rho2 = rho20 + psi2*p;

    // tmp<fvScalarMatrix> pEqnComp1;
    // tmp<fvScalarMatrix> pEqnComp2;

    // //if (transonic)
    // //{
    // //}
    // //else
    // {
    //     surfaceScalarField phid1("phid1", fvc::interpolate(psi1)*phi1);
    //     surfaceScalarField phid2("phid2", fvc::interpolate(psi2)*phi2);

    //     pEqnComp1 =
    //         fvc::ddt(rho1) + psi1*correction(fvm::ddt(p))
    //       + fvc::div(phid1, p)
    //       - fvc::Sp(fvc::div(phid1), p);

    //     pEqnComp2 =
    //         fvc::ddt(rho2) + psi2*correction(fvm::ddt(p))
    //       + fvc::div(phid2, p)
    //       - fvc::Sp(fvc::div(phid2), p);
    // }

    PtrList<surfaceScalarField> alphafs(fluid.phases().size());
    PtrList<volVectorField> HbyAs(fluid.phases().size());
    PtrList<surfaceScalarField> phiHbyAs(fluid.phases().size());
    PtrList<volScalarField> rAUs(fluid.phases().size());
    PtrList<surfaceScalarField> rAlphaAUfs(fluid.phases().size());

    phasei = 0;
    forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
    {
        phaseModel& phase = iter();

        mrfZones.makeAbsolute(phase.phi().oldTime());
        mrfZones.makeAbsolute(phase.phi());

        HbyAs.set(phasei, new volVectorField(phase.U()));
        phiHbyAs.set(phasei, new surfaceScalarField(1.0*phase.phi()));

        phasei++;
    }

    surfaceScalarField phiHbyA
    (
        IOobject
        (
            "phiHbyA",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("phiHbyA", dimArea*dimVelocity, 0)
    );

    phasei = 0;
    forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
    {
        phaseModel& phase = iter();
        const volScalarField& alpha = phase;

        alphafs.set(phasei, fvc::interpolate(alpha).ptr());
        alphafs[phasei].rename("hmm" + alpha.name());

        volScalarField dragCoeffi
        (
            IOobject
            (
                "dragCoeffi",
                runTime.timeName(),
                mesh
            ),
            fluid.dragCoeff(phase, dragCoeffs())/phase.rho(),
            zeroGradientFvPatchScalarField::typeName
        );
        dragCoeffi.correctBoundaryConditions();

        rAUs.set(phasei, (1.0/(UEqns[phasei].A() + dragCoeffi)).ptr());
        rAlphaAUfs.set
        (
            phasei,
            (
                alphafs[phasei]
               /fvc::interpolate(UEqns[phasei].A() + dragCoeffi)
            ).ptr()
        );

        HbyAs[phasei] = rAUs[phasei]*UEqns[phasei].H();

        phiHbyAs[phasei] =
        (
            (fvc::interpolate(HbyAs[phasei]) & mesh.Sf())
          + rAlphaAUfs[phasei]*fvc::ddtCorr(phase.U(), phase.phi())
        );
        mrfZones.makeRelative(phiHbyAs[phasei]);
        mrfZones.makeRelative(phase.phi().oldTime());
        mrfZones.makeRelative(phase.phi());

        phiHbyAs[phasei] +=
            rAlphaAUfs[phasei]
           *(
               fluid.surfaceTension(phase)*mesh.magSf()/phase.rho()
             + (g & mesh.Sf())
            );

        multiphaseSystem::dragModelTable::const_iterator dmIter =
            fluid.dragModels().begin();
        multiphaseSystem::dragCoeffFields::const_iterator dcIter =
            dragCoeffs().begin();
        for
        (
            ;
            dmIter != fluid.dragModels().end() && dcIter != dragCoeffs().end();
            ++dmIter, ++dcIter
        )
        {
            const phaseModel *phase2Ptr = NULL;

            if (&phase == &dmIter()->phase1())
            {
                phase2Ptr = &dmIter()->phase2();
            }
            else if (&phase == &dmIter()->phase2())
            {
                phase2Ptr = &dmIter()->phase1();
            }
            else
            {
                continue;
            }

            phiHbyAs[phasei] +=
                fvc::interpolate((*dcIter())/phase.rho())
               /fvc::interpolate(UEqns[phasei].A() + dragCoeffi)
               *phase2Ptr->phi();

            HbyAs[phasei] +=
                (1.0/phase.rho())*rAUs[phasei]*(*dcIter())
               *phase2Ptr->U();

            // Alternative flux-pressure consistent drag
            // but not momentum conservative
            //
            // HbyAs[phasei] += fvc::reconstruct
            // (
            //     fvc::interpolate
            //     (
            //         (1.0/phase.rho())*rAUs[phasei]*(*dcIter())
            //     )*phase2Ptr->phi()
            // );
        }

        phiHbyA += alphafs[phasei]*phiHbyAs[phasei];

        phasei++;
    }

    surfaceScalarField rAUf
    (
        IOobject
        (
            "rAUf",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("rAUf", dimensionSet(-1, 3, 1, 0, 0), 0)
    );

    phasei = 0;
    forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
    {
        phaseModel& phase = iter();
        rAUf += mag(alphafs[phasei]*rAlphaAUfs[phasei])/phase.rho();

        phasei++;
    }

    // Update the fixedFluxPressure BCs to ensure flux consistency
    {
        surfaceScalarField::GeometricBoundaryField phib(phi.boundaryField());
        phib = 0;
        phasei = 0;
        forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
        {
            phaseModel& phase = iter();

            phib +=
                alphafs[phasei].boundaryField()
               *(mesh.Sf().boundaryField() & phase.U().boundaryField());

            phasei++;
        }

        setSnGrad<fixedFluxPressureFvPatchScalarField>
        (
            p.boundaryField(),
            (
                phiHbyA.boundaryField() - mrfZones.relative(phib)
            )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
        );
    }

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p)
        );

        pEqnIncomp.setReference(pRefCell, pRefValue);

        solve
        (
            // (
            //     (alpha1/rho1)*pEqnComp1()
            //   + (alpha2/rho2)*pEqnComp2()
            // ) +
            pEqnIncomp,
            mesh.solver(p.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            surfaceScalarField mSfGradp("mSfGradp", pEqnIncomp.flux()/rAUf);

            phasei = 0;
            phi = dimensionedScalar("phi", phi.dimensions(), 0);
            forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
            {
                phaseModel& phase = iter();

                phase.phi() =
                    phiHbyAs[phasei]
                  + rAlphaAUfs[phasei]*mSfGradp/phase.rho();
                phi +=
                    alphafs[phasei]*phiHbyAs[phasei]
                  + mag(alphafs[phasei]*rAlphaAUfs[phasei])
                   *mSfGradp/phase.rho();

                phasei++;
            }

            // dgdt =

            // (
            //     pos(alpha2)*(pEqnComp2 & p)/rho2
            //   - pos(alpha1)*(pEqnComp1 & p)/rho1
            // );

            p.relax();
            mSfGradp = pEqnIncomp.flux()/rAUf;

            U = dimensionedVector("U", dimVelocity, vector::zero);

            phasei = 0;
            forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
            {
                phaseModel& phase = iter();
                const volScalarField& alpha = phase;

                phase.U() =
                    HbyAs[phasei]
                  + fvc::reconstruct
                    (
                        rAlphaAUfs[phasei]*(g & mesh.Sf())
                      + rAlphaAUfs[phasei]*mSfGradp/phase.rho()
                    );

                //phase.U() = fvc::reconstruct(phase.phi());
                phase.U().correctBoundaryConditions();

                U += alpha*phase.U();

                phasei++;
            }
        }
    }

    //p = max(p, pMin);

    #include "continuityErrs.H"

    // rho1 = rho10 + psi1*p;
    // rho2 = rho20 + psi2*p;

    // Dp1Dt = fvc::DDt(phi1, p);
    // Dp2Dt = fvc::DDt(phi2, p);
}
