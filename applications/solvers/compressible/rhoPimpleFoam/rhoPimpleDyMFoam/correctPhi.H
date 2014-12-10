if (mesh.changing())
{
    forAll(U.boundaryField(), patchi)
    {
        if (U.boundaryField()[patchi].fixesValue())
        {
            U.boundaryField()[patchi].initEvaluate();
        }
    }

    forAll(U.boundaryField(), patchi)
    {
        if (U.boundaryField()[patchi].fixesValue())
        {
            U.boundaryField()[patchi].evaluate();

            phi.boundaryField()[patchi] =
                rho.boundaryField()[patchi]
               *(
                   U.boundaryField()[patchi]
                 & mesh.Sf().boundaryField()[patchi]
               );
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
        dimensionedScalar("pcorr", p.dimensions(), 0.0),
        pcorrTypes
    );

    dimensionedScalar rAUf("rAUf", dimTime, 1.0);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::ddt(psi, pcorr)
          + fvc::div(phi)
          - fvm::laplacian(rAUf, pcorr)
         ==
            divrhoU
        );

        pcorrEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi += pcorrEqn.flux();
        }
    }
}
