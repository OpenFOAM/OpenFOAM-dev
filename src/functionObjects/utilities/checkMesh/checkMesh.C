/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "checkMesh.H"
#include "fvMesh.H"
#include "meshCheck.H"
#include "vtkSetWriter.H"
#include "vtkSurfaceWriter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(checkMesh, 0);
    addToRunTimeSelectionTable(functionObject, checkMesh, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::checkMesh::checkMesh
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    stopAt_(Time::stopAtControl::endTime)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::checkMesh::~checkMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::checkMesh::read(const dictionary& dict)
{
    noTopology_ =  dict.lookupOrDefault("noTopology", false);
    allGeometry_ = dict.lookupOrDefault("allGeometry", false);
    allTopology_ = dict.lookupOrDefault("allTopology", false);

    writeSurfaces_ = dict.lookupOrDefault("writeSurfaces", false);
    if (writeSurfaces_)
    {
        surfaceFormat_ = dict.lookupOrDefault<word>
        (
            "surfaceFormat",
            vtkSurfaceWriter::typeName
        );
    }

    writeSets_ = dict.lookupOrDefault("writeSets", false);
    if (writeSets_)
    {
        setFormat_ = dict.lookupOrDefault<word>
        (
            "setFormat",
            vtkSetWriter::typeName
        );
    }

    nonOrthThreshold_ = dict.lookupOrDefault("nonOrthThreshold", 70.0);
    skewThreshold_ = dict.lookupOrDefault("skewThreshold", 4.0);

    stopAt_ = Time::stopAtControlNames
    [
        dict.lookupOrDefault<word>
        (
            "stopAt",
            Time::stopAtControlNames[Time::stopAtControl::endTime]
        )
    ];

    return functionObject::read(dict);
}


bool Foam::functionObjects::checkMesh::execute()
{
    if (mesh_.changing())
    {
        autoPtr<surfaceWriter> surfWriter;
        if (writeSurfaces_)
        {
            surfWriter = surfaceWriter::New
            (
                surfaceFormat_,
                mesh_.time().writeFormat(),
                mesh_.time().writeCompression()
            );
        }

        autoPtr<setWriter> setWriter;
        if (writeSets_)
        {
            setWriter = setWriter::New
            (
                setFormat_,
                mesh_.time().writeFormat(),
                mesh_.time().writeCompression()
            );
        }

        label nFailedChecks = 0;

        if (!noTopology_)
        {
            nFailedChecks += meshCheck::checkTopology
            (
                mesh_,
                allTopology_,
                surfWriter,
                setWriter
            );
        }

        nFailedChecks += meshCheck::checkGeometry
        (
            mesh_,
            allGeometry_,
            nonOrthThreshold_,
            skewThreshold_,
            surfWriter,
            setWriter
        );

        if (nFailedChecks == 0)
        {
            Info<< "\n    Mesh OK.\n" << endl;
        }
        else
        {
            Info<< "\n    Failed " << nFailedChecks << " mesh checks.\n";

            if (stopAt_ != Time::stopAtControl::endTime)
            {
                Info<< "    Stopping at " << Time::stopAtControlNames[stopAt_]
                    << endl;

                time_.stopAt(stopAt_);
            }

            Info<< endl;
        }

        return nFailedChecks == 0;
    }
    else
    {
        return true;
    }
}


bool Foam::functionObjects::checkMesh::write()
{
    return true;
}


// ************************************************************************* //
