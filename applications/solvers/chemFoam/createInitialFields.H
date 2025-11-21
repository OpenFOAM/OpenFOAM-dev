// write initial fields - not registered since will be re-read by
// thermo package

Info<< "Reading initial conditions.\n" << endl;
IOdictionary initialConditions
(
    IOobject
    (
        "initialConditions",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const word constProp(initialConditions.lookup("constantProperty"));

if
(
    (constProp != "pressure")
 && (constProp != "volume")
 && (constProp != "temperature")
)
{
    FatalError
        << "in initialConditions, unknown constantProperty type "
        << constProp << nl
        << "    Valid types are: pressure volume temperature."
        << exit(FatalError);
}

const scalar p0 = initialConditions.lookup<scalar>("p");
const scalar T0 = initialConditions.lookup<scalar>("T");

Info<< "Creating base fields for time " << runTime.name() << endl;
{
    const word fractionBasis(initialConditions.lookup("fractionBasis"));
    if ((fractionBasis != "mass") && (fractionBasis != "mole"))
    {
        FatalError << "in initialConditions, unknown fractionBasis type " << nl
            << "Valid types are: mass or mole."
            << fractionBasis << exit(FatalError);
    }

    const dictionary specieFractions(initialConditions.subDict("fractions"));

    const label nSpecie = specieFractions.size();
    wordList specieNames(nSpecie);
    scalarList Y0(nSpecie);

    if (fractionBasis == "mole")
    {
        Info<< "Reading physicalProperties\n" << endl;
        // To obtain the specie molecular weights
        IOdictionary physicalProperties
        (
            IOobject
            (
                "physicalProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );

        scalarList W(nSpecie);
        scalarList X0(nSpecie);

        label i = 0;
        forAllConstIter(dictionary, specieFractions, iter)
        {
            specieNames[i] = iter().keyword();
            W[i] = readScalar
            (
                physicalProperties.subDict(specieNames[i]).subDict("specie")
               .lookup("molWeight")
            );
            X0[i++] = readScalar(iter().stream());
        }

        scalar mw = 0.0;
        const scalar mTot = sum(X0);
        forAll(X0, i)
        {
            X0[i] /= mTot;
            mw += W[i]*X0[i];
        }

        forAll(Y0, i)
        {
            Y0[i] = X0[i]*W[i]/mw;
        }
    }
    else  // mass fraction
    {
        label i = 0;
        forAllConstIter(dictionary, specieFractions, iter)
        {
            specieNames[i] = iter().keyword();
            Y0[i] = readScalar(iter().stream());
        }

        const scalar mTot = sum(Y0);
        forAll(Y0, i)
        {
            Y0[i] /= mTot;
        }
    }

    // Create and write the initial fields

    forAll(specieNames, i)
    {
        volScalarField
        (
            IOobject
            (
                specieNames[i],
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimless, Y0[i])
        ).write();
    }

    volScalarField
    (
        IOobject
        (
            "Ydefault",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(dimless, 0)
    ).write();

    volScalarField
    (
        IOobject
        (
            "p",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(dimPressure, p0)
    ).write();

    volScalarField
    (
        IOobject
        (
            "T",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(dimTemperature, T0)
    ).write();
}
