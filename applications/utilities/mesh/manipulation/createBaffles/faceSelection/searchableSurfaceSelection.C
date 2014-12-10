/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "searchableSurfaceSelection.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"
#include "searchableSurface.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace faceSelections
{
    defineTypeNameAndDebug(searchableSurfaceSelection, 0);
    addToRunTimeSelectionTable
    (
        faceSelection,
        searchableSurfaceSelection,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceSelections::searchableSurfaceSelection::searchableSurfaceSelection
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    faceSelection(name, mesh, dict),
    surfacePtr_
    (
        searchableSurface::New
        (
            word(dict.lookup("surface")),
            IOobject
            (
                dict.lookupOrDefault("name", mesh.objectRegistry::db().name()),
                mesh.time().constant(),
                "triSurface",
                mesh.objectRegistry::db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceSelections::searchableSurfaceSelection::~searchableSurfaceSelection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceSelections::searchableSurfaceSelection::select
(
    const label zoneID,
    labelList& faceToZoneID,
    boolList& faceToFlip
) const
{
    // Get cell-cell centre vectors

    pointField start(mesh_.nFaces());
    pointField end(mesh_.nFaces());

    // Internal faces
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        start[faceI] = mesh_.cellCentres()[mesh_.faceOwner()[faceI]];
        end[faceI] = mesh_.cellCentres()[mesh_.faceNeighbour()[faceI]];
    }

    // Boundary faces
    vectorField neighbourCellCentres;
    syncTools::swapBoundaryCellPositions
    (
        mesh_,
        mesh_.cellCentres(),
        neighbourCellCentres
    );

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start()+i;
                start[faceI] = mesh_.cellCentres()[mesh_.faceOwner()[faceI]];
                end[faceI] = neighbourCellCentres[faceI-mesh_.nInternalFaces()];
            }
        }
        else
        {
            forAll(pp, i)
            {
                label faceI = pp.start()+i;
                start[faceI] = mesh_.cellCentres()[mesh_.faceOwner()[faceI]];
                end[faceI] = mesh_.faceCentres()[faceI];
            }
        }
    }

    List<pointIndexHit> hits;
    surfacePtr_().findLine(start, end, hits);
    pointField normals;
    surfacePtr_().getNormal(hits, normals);

    //- Note: do not select boundary faces.

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        if (hits[faceI].hit())
        {
            faceToZoneID[faceI] = zoneID;
            vector d = end[faceI]-start[faceI];
            faceToFlip[faceI] = ((normals[faceI] & d) < 0);
        }
    }
    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start()+i;
                if (hits[faceI].hit())
                {
                    faceToZoneID[faceI] = zoneID;
                    vector d = end[faceI]-start[faceI];
                    faceToFlip[faceI] = ((normals[faceI] & d) < 0);
                }
            }
        }
    }

    faceSelection::select(zoneID, faceToZoneID, faceToFlip);
}


// ************************************************************************* //
