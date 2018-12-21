// Fields for data gathering

List<scalarField> allSpeciesN_RU
(
    molecules.potential().nIds(),
    scalarField(mesh.nCells(), 0.0)
);

List<scalarField> allSpeciesM_RU
(
    molecules.potential().nIds(),
    scalarField(mesh.nCells(), 0.0)
);

List<vectorField> allSpeciesVelocitySum_RU
(
    molecules.potential().nIds(),
    vectorField(mesh.nCells(), Zero)
);

List<scalarField> allSpeciesVelocityMagSquaredSum_RU
(
    molecules.potential().nIds(),
    scalarField(mesh.nCells(), 0.0)
);

// Geometric Fields for IO

Info << nl << "Creating fields." << endl;

/*---------------------------------------------------------------------------*\
    Number density
\*---------------------------------------------------------------------------*/

PtrList<volScalarField> allSpeciesRhoN
(
    molecules.potential().nIds()
);

forAll(allSpeciesRhoN, rN)
{
    Info<< "    Creating number density field for "
        << molecules.potential().idList()[rN] << endl;

    allSpeciesRhoN.set
    (
        rN,
        new volScalarField
        (
            IOobject
            (
                "rhoN_" + molecules.potential().idList()[rN],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimless/dimVolume,
            "zeroGradient"
        )
    );
    allSpeciesRhoN[rN].primitiveFieldRef() = 0;
    allSpeciesRhoN[rN].correctBoundaryConditions();
}

Info << "    Creating total number density field" << endl;

volScalarField totalRhoN
(
    IOobject
    (
        "rhoN_total",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimless/dimVolume,
    "zeroGradient"
);
totalRhoN.primitiveFieldRef() = 0;
totalRhoN.correctBoundaryConditions();

/*---------------------------------------------------------------------------*\
    Mass density
\*---------------------------------------------------------------------------*/

PtrList<volScalarField> allSpeciesRhoM
(
    molecules.potential().nIds()
);

forAll(allSpeciesRhoM, rM)
{
    Info<< "    Creating mass density field for "
        << molecules.potential().idList()[rM] << endl;

    allSpeciesRhoM.set
    (
        rM,
        new volScalarField
        (
            IOobject
            (
                "rhoM_" + molecules.potential().idList()[rM],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimMass/dimVolume,
            "zeroGradient"
        )
    );
    allSpeciesRhoM[rM].primitiveFieldRef() = 0;
    allSpeciesRhoM[rM].correctBoundaryConditions();
}

Info << "    Creating total mass density field" << endl;

volScalarField totalRhoM
(
    IOobject
    (
        "rhoM_total",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimMass/dimVolume,
    "zeroGradient"
);
totalRhoM.primitiveFieldRef() = 0;
totalRhoM.correctBoundaryConditions();

/*---------------------------------------------------------------------------*\
    Bulk velocity
\*---------------------------------------------------------------------------*/

PtrList<volVectorField> allSpeciesVelocity
(
    molecules.potential().nIds()
);

forAll(allSpeciesVelocity, v)
{
    Info<< "    Creating velocity field for "
        << molecules.potential().idList()[v] << endl;

    allSpeciesVelocity.set
    (
        v,
        new volVectorField
        (
            IOobject
            (
                "velocity_" + molecules.potential().idList()[v],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimVelocity,
            "zeroGradient"
        )
    );
    allSpeciesVelocity[v].primitiveFieldRef() = Zero;
    allSpeciesVelocity[v].correctBoundaryConditions();
}

Info << "    Creating total velocity field" << endl;

// volVectorField totalVelocity
// (
//     IOobject
//     (
//         "velocity_total",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimVelocity,
//     "zeroGradient"
// );
// totalVelocity.primitiveFieldRef() = Zero;
// totalVelocity.correctBoundaryConditions();


volVectorField totalVelocity
(
    IOobject
    (
        "velocity_total",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimVelocity, Zero)
);

/*---------------------------------------------------------------------------*\
    Kinetic temperature
\*---------------------------------------------------------------------------*/

PtrList<volScalarField> allSpeciesTemperature
(
    molecules.potential().nIds()
);

forAll(allSpeciesTemperature, t)
{
    Info<< "    Creating temperature field for "
        << molecules.potential().idList()[t] << endl;

    allSpeciesTemperature.set
    (
        t,
        new volScalarField
        (
            IOobject
            (
                "temperature_" + molecules.potential().idList()[t],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimTemperature,
            "zeroGradient"
        )
    );
    allSpeciesTemperature[t].primitiveFieldRef() = 0;
    allSpeciesTemperature[t].correctBoundaryConditions();
}

Info << "    Creating total temperature field" << endl;

volScalarField totalTemperature
(
    IOobject
    (
        "temperature_total",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimTemperature,
    "zeroGradient"
);
totalTemperature.primitiveFieldRef() = 0;
totalTemperature.correctBoundaryConditions();

/*---------------------------------------------------------------------------*\
    Mean kinetic energy
\*---------------------------------------------------------------------------*/


PtrList<volScalarField> allSpeciesMeanKE
(
    molecules.potential().nIds()
);

forAll(allSpeciesMeanKE, mKE)
{
    Info<< "    Creating mean kinetic energy field for "
        << molecules.potential().idList()[mKE] << endl;

    allSpeciesMeanKE.set
    (
        mKE,
        new volScalarField
        (
            IOobject
            (
                "meanKE_" + molecules.potential().idList()[mKE],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionSet(1, 2, -2, 0, 0, 0, 0),
            "zeroGradient"
        )
    );
    allSpeciesMeanKE[mKE].primitiveFieldRef() = 0;
    allSpeciesMeanKE[mKE].correctBoundaryConditions();
}

Info << "    Creating total mean kinetic energy field" << endl;

volScalarField totalMeanKE
(
    IOobject
    (
        "meanKE_total",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1, 2, -2, 0, 0, 0, 0),
    "zeroGradient"
);
totalMeanKE.primitiveFieldRef() = 0;
totalMeanKE.correctBoundaryConditions();
