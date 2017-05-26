/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "levelSet.H"
#include "timeSelector.H"
#include "waveAlphaFvPatchScalarField.H"
#include "waveVelocityFvPatchVectorField.H"
#include "waveSuperposition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void addWaves
(
    const waveSuperposition& waves,
    const bool liquid,
    volScalarField& alpha,
    volVectorField& U
)
{
    const scalar t = alpha.time().value();;
    const fvMesh& mesh = alpha.mesh();
    const pointMesh& pMesh = pointMesh::New(mesh);

    // Height fields
    const scalarField heightC(waves.height(t, mesh.cellCentres()));
    const scalarField heightP(waves.height(t, mesh.points()));

    // Velocity fields
    const DimensionedField<vector, volMesh>
        UGasC
        (
            IOobject("UGasC", mesh.time().timeName(), mesh),
            mesh,
            dimVelocity,
            waves.UGas(t, mesh.cellCentres())
        );
    const DimensionedField<vector, pointMesh>
        UGasP
        (
            IOobject("UGasP", mesh.time().timeName(), mesh),
            pMesh,
            dimVelocity,
            waves.UGas(t, mesh.points())
        );
    const DimensionedField<vector, volMesh>
        ULiquidC
        (
            IOobject("ULiquidC", mesh.time().timeName(), mesh),
            mesh,
            dimVelocity,
            waves.ULiquid(t, mesh.cellCentres())
        );
    const DimensionedField<vector, pointMesh>
        ULiquidP
        (
            IOobject("ULiquidP", mesh.time().timeName(), mesh),
            pMesh,
            dimVelocity,
            waves.ULiquid(t, mesh.points())
        );

    // Convert from the level set to volume-averaged fields and sum up
    alpha.ref() += levelSetFraction(mesh, heightC, heightP, !liquid);
    U.ref() +=
        levelSetAverage
        (
            mesh,
            heightC,
            heightP,
            UGasC,
            UGasP,
            ULiquidC,
            ULiquidP
        );

    // Now set the boundary fields
    forAll(alpha.boundaryField(), patchi)
    {
        fvPatchScalarField& alphap = alpha.boundaryFieldRef()[patchi];
        fvPatchVectorField& Up = U.boundaryFieldRef()[patchi];

        const fvPatch& patch = alphap.patch();

        // Height fields
        const scalarField heightF(waves.height(t, patch.Cf()));
        const scalarField heightP(waves.height(t, patch.patch().localPoints()));

        // Velocity fields
        const vectorField UGasC(waves.UGas(t, mesh.cellCentres()));
        const vectorField UGasP(waves.UGas(t, mesh.points()));
        const vectorField ULiquidC(waves.ULiquid(t, mesh.cellCentres()));
        const vectorField ULiquidP(waves.ULiquid(t, mesh.points()));

        alphap == alphap + levelSetFraction(patch, heightF, heightP, !liquid);
        Up == Up
          + levelSetAverage
            (
                patch,
                heightC,
                heightP,
                UGasC,
                UGasP,
                ULiquidC,
                ULiquidP
            );
    }
}


int main(int argc, char *argv[])
{
    timeSelector::addOptions(false, false);

    Foam::argList::addOption
    (
        "U",
        "name",
        "name of the velocity field, default is \"U\""
    );

    Foam::argList::addOption
    (
        "alpha",
        "name",
        "name of the volume fraction field, default is \"alpha\""
    );

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::selectIfPresent(runTime, args);

    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << nl << endl;

        mesh.readUpdate();

        // Read the phase fraction and velocity fields
        volScalarField alpha
        (
            IOobject
            (
                args.optionFound("alpha") ? args["alpha"] : "alpha",
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
                args.optionFound("U") ? args["U"] : "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        // Zero the fields
        alpha = dimensionedScalar("0", alpha.dimensions(), 0);
        U = dimensionedVector("0", U.dimensions(), vector::zero);
        forAll(alpha.boundaryField(), patchi)
        {
            alpha.boundaryFieldRef()[patchi] == 0;
            U.boundaryFieldRef()[patchi] == vector::zero;
        }

        // Loop the patches, looking for wave conditions
        forAll(alpha.boundaryField(), patchi)
        {
            const fvPatchScalarField& alphap = alpha.boundaryField()[patchi];
            const fvPatchVectorField& Up = U.boundaryField()[patchi];

            const bool isWave = isA<waveAlphaFvPatchScalarField>(alphap);

            if (isA<waveVelocityFvPatchVectorField>(Up) != isWave)
            {
                FatalErrorInFunction
                    << "The alpha condition on patch " << Up.patch().name()
                    << " is " << alphap.type() << " and the velocity condition"
                    << " is " << Up.type() << ". Wave boundary conditions must"
                    << " be set in pairs. If the alpha condition is "
                    << waveAlphaFvPatchScalarField::typeName
                    << " then the velocity condition must be "
                    << waveVelocityFvPatchVectorField::typeName
                    << " and vice-versa." << exit(FatalError);
            }

            if (isWave)
            {
                Info<< "Adding waves from patch " << Up.patch().name() << endl;
                addWaves
                (
                    refCast<const waveVelocityFvPatchVectorField>(Up).waves(),
                    refCast<const waveAlphaFvPatchScalarField>(alphap).liquid(),
                    alpha,
                    U
                );
            }
        }

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
