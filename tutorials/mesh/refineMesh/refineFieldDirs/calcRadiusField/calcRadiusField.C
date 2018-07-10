/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Write the volScalarField "radiusFieldXY" that has the distance to the
    origin over X,Y.

    And also write the direction fields based on the option "-calcDirections".
    The resulting fields are:
      - radialDirection
      - angularDirection
      - heightDirection

    Derived from:
      $FOAM_UTILITIES/postProcessing/miscellaneous/writeCellCentres

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
    #include "addRegionOption.H"

    argList::addBoolOption
    (
        "calcDirections",
        "calculate the direction fields as well"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);
    const bool calcDirections = args.optionFound("calcDirections");

    #include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        // Check for new mesh
        mesh.readUpdate();

        Info<< "Writing radius field over X,Y in "
            <<  runTime.timeName() << endl;

        volScalarField radiusFieldXY
        (
            IOobject
            (
                "radiusFieldXY",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sqrt
            (
                mesh.C().component(0)*mesh.C().component(0)
              + mesh.C().component(1)*mesh.C().component(1)
            )
        );
        radiusFieldXY.write();


        if(calcDirections)
        {

            vectorIOField radialDirection
            (
                IOobject
                (
                    "radialDirection",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh.C()/magSqr(mesh.C())
            );
            radialDirection.replace(vector::Z, scalar(0));
            radialDirection /= sqrt(magSqr(radialDirection));
            radialDirection.write();


            const tensor transform2Tangencial
            (
                0, -1, 0,
                1,  0, 0,
                0,  0, 1
            );
            vectorIOField angularDirection
            (
                IOobject
                (
                    "angularDirection",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                transform2Tangencial & mesh.C()
            );
            angularDirection.replace(vector::Z, scalar(0));
            angularDirection /= sqrt(magSqr(angularDirection));
            angularDirection.write();

            vectorIOField heightDirection
            (
                IOobject
                (
                    "heightDirection",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                radialDirection ^ angularDirection
            );
            heightDirection.write();
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
