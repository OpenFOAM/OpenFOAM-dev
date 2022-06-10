/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "fvMeshTopoChangersRaw.H"
#include "polyTopoChangeMap.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{
    defineTypeNameAndDebug(raw, 0);
    addToRunTimeSelectionTable(fvMeshTopoChanger, raw, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::raw::raw(fvMesh& mesh, const dictionary& dict)
:
    fvMeshTopoChanger(mesh),
    topoChanger_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::raw::~raw()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshTopoChangers::raw::update()
{
    // Do mesh changes (use inflation - put new points in topoChangeMap)
    Info<< "raw : Checking for topology changes..."
        << endl;

    // Do any topology changes. Sets topoChanged (through polyTopoChange)
    autoPtr<polyTopoChangeMap> topoChangeMap = topoChanger_.changeMesh(true);

    bool hasChanged = topoChangeMap.valid();

    if (hasChanged)
    {
        Info<< "raw : Done topology changes..."
            << endl;

        // Temporary: fix fields on patch faces created out of nothing
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Two situations:
        // - internal faces inflated out of nothing
        // - patch faces created out of previously internal faces

        // Is face mapped in any way?
        PackedBoolList mappedFace(mesh().nFaces());

        const label nOldInternal = topoChangeMap().oldPatchStarts()[0];

        const labelList& faceMap = topoChangeMap().faceMap();
        for (label facei = 0; facei < mesh().nInternalFaces(); facei++)
        {
            if (faceMap[facei] >= 0)
            {
                mappedFace[facei] = 1;
            }
        }
        for
        (
            label facei = mesh().nInternalFaces();
            facei < mesh().nFaces();
            facei++
        )
        {
            if (faceMap[facei] >= 0 && faceMap[facei] >= nOldInternal)
            {
                mappedFace[facei] = 1;
            }
        }

        const List<objectMap>& fromFaces = topoChangeMap().facesFromFacesMap();

        forAll(fromFaces, i)
        {
            mappedFace[fromFaces[i].index()] = 1;
        }

        const List<objectMap>& fromEdges = topoChangeMap().facesFromEdgesMap();

        forAll(fromEdges, i)
        {
            mappedFace[fromEdges[i].index()] = 1;
        }

        const List<objectMap>& fromPts = topoChangeMap().facesFromPointsMap();

        forAll(fromPts, i)
        {
            mappedFace[fromPts[i].index()] = 1;
        }

        // Set unmapped faces to zero
        Info<< "fvMeshTopoChangers::raw :"
            << " zeroing unmapped boundary values." << endl;
        zeroUnmappedValues<scalar, fvPatchField, volMesh>(mappedFace);
        zeroUnmappedValues<vector, fvPatchField, volMesh>(mappedFace);
        zeroUnmappedValues<sphericalTensor, fvPatchField, volMesh>(mappedFace);
        zeroUnmappedValues<symmTensor, fvPatchField, volMesh>(mappedFace);
        zeroUnmappedValues<tensor, fvPatchField, volMesh>(mappedFace);

        if (topoChangeMap().hasMotionPoints())
        {
            pointField newPoints = topoChangeMap().preMotionPoints();

            // Give the meshModifiers opportunity to modify points
            Info<< "fvMeshTopoChangers::raw :"
                << " calling modifyMotionPoints." << endl;
            topoChanger_.modifyMotionPoints(newPoints);

            // Actually move points
            Info<< "fvMeshTopoChangers::raw :"
                << " calling movePoints." << endl;

            mesh().movePoints(newPoints);
        }
    }
    else
    {
        // Pout<< "fvMeshTopoChangers::raw :"
        //    << " no topology changes..." << endl;
    }

    return hasChanged;
}


void Foam::fvMeshTopoChangers::raw::topoChange(const polyTopoChangeMap& map)
{}


void Foam::fvMeshTopoChangers::raw::mapMesh(const polyMeshMap& map)
{}


void Foam::fvMeshTopoChangers::raw::distribute
(
    const polyDistributionMap& map
)
{}


// ************************************************************************* //
