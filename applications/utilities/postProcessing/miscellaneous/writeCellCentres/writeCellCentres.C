/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Description
    Write the three components of the cell centres as volScalarFields so
    they can be used in postprocessing in thresholding.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "vectorIOField.H"
#include "volFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    timeSelector::addOptions();
#   include "addRegionOption.H"

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        // Check for new mesh
        mesh.readUpdate();

        //volVectorField cc
        //(
        //    IOobject
        //    (
        //        "cellCentres",
        //        runTime.timeName(),
        //        mesh,
        //        IOobject::NO_READ,
        //        IOobject::AUTO_WRITE
        //    ),
        //    1.0*mesh.C()
        //);
        //
        //Info<< "Writing cellCentre positions to " << cc.name() << " in "
        //     << runTime.timeName() << endl;
        //cc.write();

        Info<< "Writing components of cellCentre positions to volScalarFields"
            << " ccx, ccy, ccz in " <<  runTime.timeName() << endl;

        for (direction i=0; i<vector::nComponents; i++)
        {
            volScalarField cci
            (
                IOobject
                (
                    "cc" + word(vector::componentNames[i]),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh.C().component(i)
            );

            cci.write();
        }


        volScalarField V
        (
            IOobject
            (
                "V",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("V", mesh.V().dimensions(), 0.0),
            calculatedFvPatchField<scalar>::typeName
        );
        V.dimensionedInternalField() = mesh.V();
        forAll(V.boundaryField(), patchI)
        {
            V.boundaryField()[patchI] =
                V.boundaryField()[patchI].patch().magSf();
        }
        Info<< "Writing cellVolumes and patch faceAreas to " << V.name()
            << " in " << runTime.timeName() << endl;
        V.write();

    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
