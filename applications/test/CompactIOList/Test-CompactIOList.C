/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    testCompactIOList

Description
    Simple demonstration and test application for the CompactIOList container

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "argList.H"
#include "Time.H"
#include "polyMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//  Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    IOstream::streamFormat format=IOstream::BINARY;
    // IOstream::streamFormat format=IOstream::ASCII;

    const label size = 20000000;

    // Old format
    // ~~~~~~~~~~

    {
        // Construct big faceList in old format
        faceIOList faces2
        (
            IOobject
            (
                "faces2",
                runTime.constant(),
                polyMesh::meshSubDir,
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            size
        );

        const face f(identity(4));

        forAll(faces2, i)
        {
            faces2[i] = f;
        }

        Info<< "Constructed faceList in = "
            << runTime.cpuTimeIncrement() << " s" << nl << endl;


        // Write binary
        faces2.writeObject
        (
            format,
            IOstream::currentVersion,
            IOstream::UNCOMPRESSED
        );

        Info<< "Written old format faceList in = "
            << runTime.cpuTimeIncrement() << " s" << nl << endl;

        // Read
        faceIOList faces3
        (
            IOobject
            (
                "faces2",
                runTime.constant(),
                polyMesh::meshSubDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        Info<< "Read old format " << faces3.size() << " faceList in = "
            << runTime.cpuTimeIncrement() << " s" << nl << endl;
    }


    // New format
    // ~~~~~~~~~~

    {
        // Construct big faceList in new format
        faceCompactIOList faces2
        (
            IOobject
            (
                "faces2",
                runTime.constant(),
                polyMesh::meshSubDir,
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            size
        );

        const face f(identity(4));

        forAll(faces2, i)
        {
            faces2[i] = f;
        }

        Info<< "Constructed new format faceList in = "
            << runTime.cpuTimeIncrement() << " s" << nl << endl;


        // Write binary
        faces2.writeObject
        (
            format,
            IOstream::currentVersion,
            IOstream::UNCOMPRESSED
        );

        Info<< "Written new format faceList in = "
            << runTime.cpuTimeIncrement() << " s" << nl << endl;

        // Read
        faceCompactIOList faces3
        (
            IOobject
            (
                "faces2",
                runTime.constant(),
                polyMesh::meshSubDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        Info<< "Read new format " << faces3.size() << " faceList in = "
            << runTime.cpuTimeIncrement() << " s" << nl << endl;
    }

    return 0;
}


// ************************************************************************* //
