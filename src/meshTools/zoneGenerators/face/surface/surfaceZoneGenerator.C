/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "surfaceZoneGenerator.H"
#include "Time.H"
#include "polyMesh.H"
#include "searchableSurface.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(surface, 0);
        addToRunTimeSelectionTable
        (
            zoneGenerator,
            surface,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::surface::surface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    surfacePtr_
    (
        searchableSurface::New
        (
            word(dict.lookup("surface")),
            IOobject
            (
                dict.lookupOrDefault
                (
                    "surfaceName",
                    mesh.objectRegistry::db().name()
                ),
                mesh.time().constant(),
                searchableSurface::geometryDir(mesh.time()),
                mesh.objectRegistry::db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict
        )
    ),
    tol_(dict.lookupOrDefault<scalar>("tol", dimless, rootSmall))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::surface::~surface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::surface::generate() const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    // Get a list of "interior" faces; i.e., either internal or coupled
    labelList interiorFaceFaces(mesh_.nFaces());
    {
        forAll(faceNeighbour, facei)
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
                    nInteriorFaces++;
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

            ownCc[interiorFacei] = cc[faceOwner[facei]];
            nbrCc[interiorFacei] =
                mesh_.isInternalFace(facei)
              ? cc[faceNeighbour[facei]]
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

            side[faceOwner[facei]] = sign ? -1 : +1;

            if (mesh_.isInternalFace(facei))
            {
                side[faceNeighbour[facei]] = sign ? +1 : -1;
            }
        }
    }
    labelList boundaryNbrSide;
    syncTools::swapBoundaryCellList(mesh_, side, boundaryNbrSide);

    // Select intersected faces and set orientation (flipMap)
    labelList faceIndices(interiorFaceFaces.size());
    boolList flipMap(interiorFaceFaces.size());

    label fi = 0;
    forAll(interiorFaceFaces, interiorFacei)
    {
        const label facei = interiorFaceFaces[interiorFacei];

        const label ownSide = side[faceOwner[facei]];
        const label nbrSide =
            mesh_.isInternalFace(facei)
          ? side[faceNeighbour[facei]]
          : boundaryNbrSide[facei - mesh_.nInternalFaces()];

        if ((ownSide == 1 && nbrSide == -1) || (ownSide == -1 && nbrSide == 1))
        {
            faceIndices[fi] = facei;
            flipMap[fi] = nbrSide > ownSide;
            fi++;
        }
    }

    faceIndices.setSize(fi);
    flipMap.setSize(fi);

    return zoneSet
    (
        new faceZone
        (
            zoneName_,
            faceIndices,
            flipMap,
            mesh_.faceZones(),
            moveUpdate_,
            true
        )
    );
}


// ************************************************************************* //
