/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    setAtmBoundaryLayer

Description
    Applies atmospheric boundary layer models to the entire domain for case
    initialisation.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "volFields.H"
#include "wallPolyPatch.H"
#include "atmBoundaryLayer.H"
#include "systemDict.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions(false, false);

    #include "addDictOption.H"
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::selectIfPresent(runTime, args);

    #include "createNamedMesh.H"

    const dictionary setAtmBoundaryLayerDict
    (
        systemDict("setAtmBoundaryLayerDict", args, mesh)
    );

    // Parse options
    const word UName = setAtmBoundaryLayerDict.lookupOrDefault<word>("U", "U");
    const word kName = setAtmBoundaryLayerDict.lookupOrDefault<word>("k", "k");
    const word epsilonName =
        setAtmBoundaryLayerDict.lookupOrDefault<word>("epsilon", "epsilon");

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        mesh.readUpdate();

        // Read the fields
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
        IOobject kIO
        (
            kName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );
        tmp<volScalarField> k
        (
            kIO.headerOk() ? new volScalarField(kIO, mesh) : nullptr
        );
        IOobject epsilonIO
        (
            epsilonName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );
        tmp<volScalarField> epsilon
        (
            epsilonIO.headerOk() ? new volScalarField(epsilonIO, mesh) : nullptr
        );

        // Set the internal fields
        const atmBoundaryLayer atm(mesh.C(), setAtmBoundaryLayerDict);
        U.primitiveFieldRef() = atm.U(mesh.C());
        if (k.valid())
        {
            k.ref().primitiveFieldRef() = atm.k(mesh.C());
        }
        if (epsilon.valid())
        {
            epsilon.ref().primitiveFieldRef() = atm.epsilon(mesh.C());
        }

        // Set all non-wall patch fields
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& patch = mesh.boundary()[patchi];

            const atmBoundaryLayer atmp(patch.Cf(), setAtmBoundaryLayerDict);

            if (!isA<wallPolyPatch>(patch))
            {
                U.boundaryFieldRef()[patchi] == atmp.U(patch.Cf());
                if (k.valid())
                {
                    k.ref().boundaryFieldRef()[patchi] ==
                        atmp.k(patch.Cf());
                }
                if (epsilon.valid())
                {
                    epsilon.ref().boundaryFieldRef()[patchi] ==
                        atmp.epsilon(patch.Cf());
                }
            }
        }

        // Output
        Info<< "Writing " << U.name() << nl << endl;
        U.write();
        if (k.valid())
        {
            Info<< "Writing " << k->name() << nl << endl;
            k->write();
        }
        if (epsilon.valid())
        {
            Info<< "Writing " << epsilon->name() << nl << endl;
            epsilon->write();
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
