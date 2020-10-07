/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    setWaves

Description
    Applies wave models to the entire domain for case initialisation using
    level sets for second-order accuracy.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "levelSet.H"
#include "pointFields.H"
#include "timeSelector.H"
#include "uniformDimensionedFields.H"
#include "volFields.H"
#include "wallPolyPatch.H"
#include "waveAlphaFvPatchScalarField.H"
#include "waveVelocityFvPatchVectorField.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions(false, false);

    #include "addDictOption.H"
    #include "addRegionOption.H"

    argList::addOption
    (
        "alpha",
        "name",
        "name of the volume fraction field, default is \"alpha\""
    );

    argList::addOption
    (
        "U",
        "name",
        "name of the velocity field, default is \"U\""
    );

    argList::addBoolOption
    (
        "gas",
        "the volume fraction field is that of the gas above the wave"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::selectIfPresent(runTime, args);

    #include "createNamedMesh.H"

    const word dictName("setWavesDict");
    #include "setSystemMeshDictionaryIO.H"
    Info<< "Reading " << dictName << "\n" << endl;
    IOdictionary setWavesDict(dictIO);

    #include "readGravitationalAcceleration.H"

    const pointMesh& pMesh = pointMesh::New(mesh);

    // Parse options
    const word alphaName = setWavesDict.lookupOrDefault<word>("alpha", "alpha");
    const word UName = setWavesDict.lookupOrDefault<word>("U", "U");
    const bool liquid = setWavesDict.lookupOrDefault<bool>("liquid", true);

    // Get the wave models
    const waveSuperposition& waves = waveSuperposition::New(mesh);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        const scalar t = runTime.value();

        Info<< "Time = " << runTime.timeName() << nl << endl;

        mesh.readUpdate();

        // Read the fields which are to be set
        volScalarField alpha
        (
            IOobject
            (
                alphaName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );
        volVectorField U
        (
            IOobject
            (
                UName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        // Create modelled fields on both cells and points
        volScalarField h
        (
            IOobject("h", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar(dimLength, 0)
        );
        pointScalarField hp
        (
            IOobject("hp", runTime.timeName(), mesh),
            pMesh,
            dimensionedScalar(dimLength, 0)
        );
        volVectorField uGas
        (
            IOobject("uGas", runTime.timeName(), mesh),
            mesh,
            dimensionedVector(dimVelocity, vector::zero)
        );
        pointVectorField uGasp
        (
            IOobject("uGasp", runTime.timeName(), mesh),
            pMesh,
            dimensionedVector(dimVelocity, vector::zero)
        );
        volVectorField uLiq
        (
            IOobject("uLiq", runTime.timeName(), mesh),
            mesh,
            dimensionedVector(dimVelocity, vector::zero)
        );
        pointVectorField uLiqp
        (
            IOobject("uLiqp", runTime.timeName(), mesh),
            pMesh,
            dimensionedVector(dimVelocity, vector::zero)
        );

        // Cell centres and points
        const pointField& ccs = mesh.cellCentres();
        const pointField& pts = mesh.points();

        // Internal field
        h.primitiveFieldRef() = waves.height(t, ccs);
        hp.primitiveFieldRef() = waves.height(t, pts);
        uGas.primitiveFieldRef() = waves.UGas(t, ccs);
        uGasp.primitiveFieldRef() = waves.UGas(t, pts);
        uLiq.primitiveFieldRef() = waves.ULiquid(t, ccs);
        uLiqp.primitiveFieldRef() = waves.ULiquid(t, pts);

        // Boundary fields
        forAll(mesh.boundary(), patchj)
        {
            const pointField& fcs = mesh.boundary()[patchj].Cf();
            h.boundaryFieldRef()[patchj] = waves.height(t, fcs);
            uGas.boundaryFieldRef()[patchj] = waves.UGas(t, fcs);
            uLiq.boundaryFieldRef()[patchj] = waves.ULiquid(t, fcs);
        }

        // Calculate the fields
        volScalarField alphaNoBCs(levelSetFraction(h, hp, !liquid));
        volVectorField UNoBCs(levelSetAverage(h, hp, uGas, uGasp, uLiq, uLiqp));

        // Set the wave and non-wall fixed-value patch fields
        forAll(mesh.boundary(), patchi)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchi];

            fvPatchScalarField& alphap = alpha.boundaryFieldRef()[patchi];
            fvPatchVectorField& Up = U.boundaryFieldRef()[patchi];
            if
            (
               !isA<wallPolyPatch>(patch)
             || isA<waveAlphaFvPatchScalarField>(alphap)
             || isA<waveVelocityFvPatchVectorField>(Up)
            )
            {
                alphap == alphaNoBCs.boundaryField()[patchi];
                Up == UNoBCs.boundaryField()[patchi];
            }
        }

        // Set the internal fields and all non-fixed value patch fields
        alpha = alphaNoBCs;
        U = UNoBCs;

        // Output
        Info<< "Writing " << alpha.name() << nl;
        alpha.write();
        Info<< "Writing " << U.name() << nl << endl;
        U.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
