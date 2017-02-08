        volScalarField::Internal Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("Sp", dgdt.dimensions(), 0.0)
        );

        volScalarField::Internal Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            divU*min(alpha1, scalar(1))
        );

        forAll(dgdt, celli)
        {
            if (dgdt[celli] > 0.0 && alpha1[celli] > 0.0)
            {
                Sp[celli] -= dgdt[celli]*alpha1[celli];
                Su[celli] += dgdt[celli]*alpha1[celli];
            }
            else if (dgdt[celli] < 0.0 && alpha1[celli] < 1.0)
            {
                Sp[celli] += dgdt[celli]*(1.0 - alpha1[celli]);
            }
        }
