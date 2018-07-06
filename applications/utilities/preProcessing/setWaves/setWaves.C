/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "fvCFD.H"
#include "levelSet.H"
#include "pointFields.H"
#include "timeSelector.H"
#include "wallDist.H"
#include "waveAlphaFvPatchScalarField.H"
#include "waveVelocityFvPatchVectorField.H"
#include "waveSuperposition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

    const pointMesh& pMesh = pointMesh::New(mesh);

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

        // Create modelled fields on both cells and points
        volScalarField h
        (
            IOobject("h", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar("0", dimLength, 0)
        );
        pointScalarField hp
        (
            IOobject("hp", runTime.timeName(), mesh),
            pMesh,
            dimensionedScalar("0", dimLength, 0)
        );
        volVectorField uGas
        (
            IOobject("uGas", runTime.timeName(), mesh),
            mesh,
            dimensionedVector("0", dimVelocity, vector::zero)
        );
        pointVectorField uGasp
        (
            IOobject("uGasp", runTime.timeName(), mesh),
            pMesh,
            dimensionedVector("0", dimLength, vector::zero)
        );
        volVectorField uLiq
        (
            IOobject("uLiq", runTime.timeName(), mesh),
            mesh,
            dimensionedVector("0", dimVelocity, vector::zero)
        );
        pointVectorField uLiqp
        (
            IOobject("uLiqp", runTime.timeName(), mesh),
            pMesh,
            dimensionedVector("0", dimLength, vector::zero)
        );

        // The number of wave patches
        label nWaves = 0;

        // Whether the alpha conditions refer to the liquid phase
        bool liquid = false;

        // Loop the patches, averaging and superimposing wave model data
        forAll(mesh.boundary(), patchi)
        {
            fvPatchScalarField& alphap = alpha.boundaryFieldRef()[patchi];
            fvPatchVectorField& Up = U.boundaryFieldRef()[patchi];

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

            if (!isWave)
            {
                continue;
            }

            Info<< "Adding waves from patch " << Up.patch().name() << endl;

            const waveSuperposition& waves =
                refCast<waveVelocityFvPatchVectorField>(Up).waves();

            const bool liquidp =
                refCast<waveAlphaFvPatchScalarField>(alphap).liquid();
            if (nWaves > 0 && liquidp != liquid)
            {
                FatalErrorInFunction
                    << "All " << waveAlphaFvPatchScalarField::typeName
                    << "patch fields must be configured for the same phase,"
                    << " i.e., the liquid switch must have the same value."
                    << exit(FatalError);
            }
            liquid = liquidp;

            const pointField& ccs = mesh.cellCentres();
            const pointField& pts = mesh.points();

            // Internal field superposition
            h.primitiveFieldRef() += waves.height(t, ccs);
            hp.primitiveFieldRef() += waves.height(t, pts);
            uGas.primitiveFieldRef() += waves.UGas(t, ccs) - waves.UMean(t);
            uGasp.primitiveFieldRef() += waves.UGas(t, pts) - waves.UMean(t);
            uLiq.primitiveFieldRef() += waves.ULiquid(t, ccs) - waves.UMean(t);
            uLiqp.primitiveFieldRef() += waves.ULiquid(t, pts) - waves.UMean(t);

            // Boundary field superposition
            forAll(mesh.boundary(), patchj)
            {
                const pointField& fcs = mesh.boundary()[patchj].Cf();
                h.boundaryFieldRef()[patchj] += waves.height(t, fcs);
                uGas.boundaryFieldRef()[patchj] +=
                    waves.UGas(t, fcs) - waves.UMean(t);
                uLiq.boundaryFieldRef()[patchj] +=
                    waves.ULiquid(t, fcs) - waves.UMean(t);
            }

            ++ nWaves;
        }

        // Create the mean velocity field
        volVectorField UMean
        (
            IOobject("UMean", runTime.timeName(), mesh),
            mesh,
            dimensionedVector("UMean", dimVelocity, Zero)
        );

        if (nWaves == 0)
        {
            // Warn and skip to the next time if there are no wave patches
            WarningInFunction
                << "No " << waveAlphaFvPatchScalarField::typeName << " or "
                << waveVelocityFvPatchVectorField::typeName << " patch fields "
                << "were found. No waves have been set." << endl;

            continue;
        }
        else if (nWaves == 1)
        {
            // Set a mean velocity equal to that on the only wave patch
            forAll(mesh.boundary(), patchi)
            {
                const fvPatchVectorField& Up = U.boundaryField()[patchi];
                if (!isA<waveVelocityFvPatchVectorField>(Up))
                {
                    continue;
                }

                const waveSuperposition& waves =
                    refCast<const waveVelocityFvPatchVectorField>(Up).waves();

                UMean ==
                    dimensionedVector("UMean", dimVelocity, waves.UMean(t));
            }
        }
        else if (nWaves > 1)
        {
            // Set the mean velocity by distance weighting from the wave patches

            // Create weighted average fields for the mean velocity
            volScalarField weight
            (
                IOobject("weight", runTime.timeName(), mesh),
                mesh,
                dimensionedScalar("0", dimless/dimLength, 0)
            );
            volVectorField weightUMean
            (
                IOobject("weightUMean", runTime.timeName(), mesh),
                mesh,
                dimensionedVector("0", dimVelocity/dimLength, vector::zero)
            );

            // Loop the patches, inverse-distance weighting the mean velocities
            forAll(mesh.boundary(), patchi)
            {
                const fvPatchVectorField& Up = U.boundaryField()[patchi];
                if (!isA<waveVelocityFvPatchVectorField>(Up))
                {
                    continue;
                }

                const waveSuperposition& waves =
                    refCast<const waveVelocityFvPatchVectorField>(Up).waves();

                const volScalarField w
                (
                    1
                   /(
                       wallDist(mesh, labelList(1, patchi)).y()
                     + dimensionedScalar("ySmall", dimLength, small)
                    )
                );

                weight += w;
                weightUMean +=
                    w*dimensionedVector("wUMean", dimVelocity, waves.UMean(t));
            }

            // Complete the average for the mean velocity
            UMean = weightUMean/weight;
        }

        // Set the fields
        alpha == levelSetFraction(h, hp, !liquid);
        U == UMean + levelSetAverage(h, hp, uGas, uGasp, uLiq, uLiqp);

        // Set the boundary fields
        forAll(mesh.boundary(), patchi)
        {
            fvPatchScalarField& alphap = alpha.boundaryFieldRef()[patchi];
            fvPatchVectorField& Up = U.boundaryFieldRef()[patchi];
            if (isA<waveAlphaFvPatchScalarField>(alphap))
            {
                alphap == refCast<waveAlphaFvPatchScalarField>(alphap).alpha();
                Up == refCast<waveVelocityFvPatchVectorField>(Up).U();
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
