if (mesh.changing())
{
    forAll(U.boundaryField(), patchI)
    {
        if (U.boundaryField()[patchI].fixesValue())
        {
            U.boundaryField()[patchI].initEvaluate();
        }
    }

    forAll(U.boundaryField(), patchI)
    {
        if (U.boundaryField()[patchI].fixesValue())
        {
            U.boundaryField()[patchI].evaluate();

            phi.boundaryField()[patchI] =
                U.boundaryField()[patchI]
              & mesh.Sf().boundaryField()[patchI];
        }
    }
}

{
    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", p_rgh.dimensions(), 0.0),
        pcorrTypes
    );

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    #ifndef divUCorr
    #define divUCorr
    #endif

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::div(phi) divUCorr
        );

        pcorrEqn.setReference(pRefCell, pRefValue);
        pcorrEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux();
        }
    }

    #undef divUCorr

    #include "continuityErrs.H"
}
