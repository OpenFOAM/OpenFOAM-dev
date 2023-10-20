/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "checkMeshQuality.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "faceSet.H"
#include "meshCheck.H"
#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::checkMeshQuality
(
    const polyMesh& mesh,
    const dictionary& dict,
    const autoPtr<surfaceWriter>& writer
)
{
    label noFailedChecks = 0;

    {
        faceSet faces(mesh, "meshQualityFaces", mesh.nFaces()/100+1);
        meshCheck::checkMesh(false, mesh, dict, faces);

        label nFaces = returnReduce(faces.size(), sumOp<label>());

        if (nFaces > 0)
        {
            noFailedChecks++;

            Info<< "  <<Writing " << nFaces
                << " faces in error to set " << faces.name() << endl;
            faces.instance() = mesh.pointsInstance();
            faces.write();
            if (writer.valid())
            {
                meshCheck::mergeAndWrite(writer(), faces);
            }
        }
    }

    return noFailedChecks;
}


// ************************************************************************* //
