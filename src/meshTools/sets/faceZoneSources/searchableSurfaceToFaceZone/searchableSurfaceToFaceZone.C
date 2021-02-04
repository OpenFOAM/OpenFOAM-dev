/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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

#include "searchableSurfaceToFaceZone.H"
#include "polyMesh.H"
#include "faceZoneSet.H"
#include "searchableSurface.H"
#include "syncTools.H"
#include "Time.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableSurfaceToFaceZone, 0);
    addToRunTimeSelectionTable
    (
        topoSetSource,
        searchableSurfaceToFaceZone,
        word
    );
}


Foam::topoSetSource::addToUsageTable Foam::searchableSurfaceToFaceZone::usage_
(
    searchableSurfaceToFaceZone::typeName,
    "\n    Usage: searchableSurfaceToFaceZone surface\n\n"
    "    Select all faces whose cell-cell centre vector intersects the surface "
    "\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::searchableSurfaceToFaceZone::combine
(
    faceZoneSet& fzSet,
    const bool add
) const
{
    // Get a list of "interior" faces; i.e., either internal or coupled
    labelList interiorFaceFaces(mesh_.nFaces());
    {
        forAll(mesh_.faceNeighbour(), facei)
        {
            interiorFaceFaces[facei] = facei;
        }
        label nInteriorFaces = mesh_.nInternalFaces();
        forAll(mesh_.boundaryMesh(), patchi)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchi];
            if (patch.coupled())
            {
                forAll(patch, patchFacei)
                {
                    const label facei = patch.start() + patchFacei;
                    interiorFaceFaces[nInteriorFaces] = facei;
                    ++ nInteriorFaces;
                }
            }
        }
        interiorFaceFaces.resize(nInteriorFaces);
    }

    // Create owner and neighbour cell centre lists for all interior faces
    pointField ownCc(interiorFaceFaces.size());
    pointField nbrCc(interiorFaceFaces.size());
    {
        const pointField& cc = mesh_.cellCentres();

        vectorField boundaryNbrCc;
        syncTools::swapBoundaryCellPositions(mesh_, cc, boundaryNbrCc);

        forAll(ownCc, interiorFacei)
        {
            const label facei = interiorFaceFaces[interiorFacei];

            ownCc[interiorFacei] = cc[mesh_.faceOwner()[facei]];
            nbrCc[interiorFacei] =
                mesh_.isInternalFace(facei)
              ? cc[mesh_.faceNeighbour()[facei]]
              : boundaryNbrCc[facei - mesh_.nInternalFaces()];
        }
    }

    // Do intersection tests on the vectors between the owner and neighbour
    // cell centres, extended by the tolerance
    List<pointIndexHit> hits;
    pointField normals;
    surfacePtr_().findLine
    (
        ownCc + tol_*(ownCc - nbrCc),
        nbrCc - tol_*(ownCc - nbrCc),
        hits
    );
    surfacePtr_().getNormal(hits, normals);

    // Create a list of labels indicating what side of the surface a cell
    // is on; -1 is below, +1 is above, and 0 is too far from the surface
    // for the sidedness to be calculable
    labelList side(mesh_.nCells(), 0);
    forAll(hits, interiorFacei)
    {
        if (hits[interiorFacei].hit())
        {
            const label facei = interiorFaceFaces[interiorFacei];

            const vector d = nbrCc[interiorFacei] - ownCc[interiorFacei];
            const bool sign = (normals[interiorFacei] & d) < 0;

            side[mesh_.faceOwner()[facei]] = sign ? -1 : +1;

            if (mesh_.isInternalFace(facei))
            {
                side[mesh_.faceNeighbour()[facei]] = sign ? +1 : -1;
            }
        }
    }
    labelList boundaryNbrSide;
    syncTools::swapBoundaryCellList(mesh_, side, boundaryNbrSide);

    // Create a face direction list indicating which direction a given face
    // intersects the surface; -1 is backward, +1 is forward, and 0
    // indicates that the face does not intersect the surface
    labelList direction(mesh_.nFaces(), 0);
    forAll(interiorFaceFaces, interiorFacei)
    {
        const label facei = interiorFaceFaces[interiorFacei];

        const label ownSide = side[mesh_.faceOwner()[facei]];
        const label nbrSide =
            mesh_.isInternalFace(facei)
          ? side[mesh_.faceNeighbour()[facei]]
          : boundaryNbrSide[facei - mesh_.nInternalFaces()];

        direction[facei] = ownSide*nbrSide < 0 ? ownSide - nbrSide : 0;
    }

    // Select intersected faces
    DynamicList<label> newAddressing;
    DynamicList<bool> newFlipMap;
    if (add)
    {
        // Start from copy
        newAddressing = DynamicList<label>(fzSet.addressing());
        newFlipMap = DynamicList<bool>(fzSet.flipMap());

        // Add everything with a direction that is not already in the set
        forAll(direction, facei)
        {
            if (direction[facei] != 0 && !fzSet.found(facei))
            {
                newAddressing.append(facei);
                newFlipMap.append(direction[facei] < 0);
            }
        }
    }
    else
    {
        // Start from empty
        newAddressing = DynamicList<label>(fzSet.addressing().size());
        newFlipMap = DynamicList<bool>(fzSet.flipMap().size());

        // Add everything from the zone set that does not have a direction
        forAll(fzSet.addressing(), i)
        {
            if (direction[fzSet.addressing()[i]] == 0)
            {
                newAddressing.append(fzSet.addressing()[i]);
                newFlipMap.append(fzSet.flipMap()[i]);
            }
        }
    }
    fzSet.addressing().transfer(newAddressing);
    fzSet.flipMap().transfer(newFlipMap);
    fzSet.updateSet();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaceToFaceZone::searchableSurfaceToFaceZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    surfacePtr_
    (
        searchableSurface::New
        (
            word(dict.lookup("surface")),
            IOobject
            (
                dict.lookupOrDefault("name", mesh.objectRegistry::db().name()),
                mesh.time().constant(),
                searchableSurface::geometryDir(mesh.time()),
                mesh.objectRegistry::db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict
        )
    ),
    tol_(dict.lookupOrDefault<scalar>("tol", rootSmall))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableSurfaceToFaceZone::~searchableSurfaceToFaceZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::searchableSurfaceToFaceZone::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (!isA<faceZoneSet>(set))
    {
        WarningInFunction
            << "Operation only allowed on a faceZoneSet." << endl;
    }
    else
    {
        faceZoneSet& fzSet = refCast<faceZoneSet>(set);

        if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
        {
            Info<< "    Adding all faces from surface "
                << surfacePtr_().name() << " ..." << endl;

            combine(fzSet, true);
        }
        else if (action == topoSetSource::DELETE)
        {
            Info<< "    Removing all faces from surface "
                << surfacePtr_().name() << " ..." << endl;

            combine(fzSet, false);
        }
    }
}


// ************************************************************************* //
