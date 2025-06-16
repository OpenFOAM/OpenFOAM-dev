/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "faceZoneGenerator.H"
#include "polyMesh.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(faceZoneGenerator, 0);
        addToRunTimeSelectionTable
        (
            zoneGenerator,
            faceZoneGenerator,
            dictionary
        );
    }
}

const Foam::NamedEnum
<
    Foam::zoneGenerators::faceZoneGenerator::cellFaces,
    4
> Foam::zoneGenerators::faceZoneGenerator::cellFacesNames
{
    "all",
    "inner",
    "outer",
    "outerInternal"
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::faceZoneGenerator::faceZoneGenerator
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    zoneGenerators_(mesh, dict),
    cellFaces_
    (
        cellFacesNames.lookupOrDefault
        (
            "cellFaces",
            dict,
            cellFaces::all
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::faceZoneGenerator::~faceZoneGenerator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::faceZoneGenerator::generate() const
{
    boolList selectedFaces(mesh_.nFaces(), false);
    boolList flipMap(mesh_.nFaces(), false);

    forAll(zoneGenerators_, i)
    {
        zoneSet zs(zoneGenerators_[i].generate());

        if (zs.pZone.valid())
        {
            // const labelList& zonePoints = zs.pZone();
            NotImplemented;
        }

        if (zs.cZone.valid())
        {
            const labelList& zoneCells = zs.cZone();

            switch (cellFaces_)
            {
                case cellFaces::all:
                {
                    forAll(zoneCells, zci)
                    {
                        const labelList& cellFaces =
                            mesh_.cells()[zoneCells[zci]];

                        forAll(cellFaces, cFacei)
                        {
                            selectedFaces[cellFaces[cFacei]] = true;
                        }
                    }
                    break;
                }

                case cellFaces::inner:
                {
                    const labelHashSet cellSet(zoneCells);

                    const label nInt = mesh_.nInternalFaces();
                    const labelList& own = mesh_.faceOwner();
                    const labelList& nei = mesh_.faceNeighbour();
                    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

                    // Check all internal faces
                    for (label facei = 0; facei < nInt; facei++)
                    {
                        if
                        (
                            cellSet.found(own[facei])
                         && cellSet.found(nei[facei])
                        )
                        {
                            selectedFaces[facei] = true;
                        }
                    }


                    // Get coupled cell status
                    boolList neiInSet(mesh_.nFaces()-nInt, false);

                    forAll(patches, patchi)
                    {
                        const polyPatch& pp = patches[patchi];

                        if (pp.coupled())
                        {
                            label facei = pp.start();
                            forAll(pp, i)
                            {
                                neiInSet[facei-nInt] =
                                    cellSet.found(own[facei]);
                                facei++;
                            }
                        }
                    }
                    syncTools::swapBoundaryFaceList(mesh_, neiInSet);

                    // Check all boundary faces
                    forAll(patches, patchi)
                    {
                        const polyPatch& pp = patches[patchi];

                        if (pp.coupled())
                        {
                            label facei = pp.start();
                            forAll(pp, i)
                            {
                                if
                                (
                                    cellSet.found(own[facei])
                                 && neiInSet[facei-nInt]
                                )
                                {
                                    selectedFaces[facei] = true;
                                }
                                facei++;
                            }
                        }
                    }
                    break;
                }

                case cellFaces::outer:
                case cellFaces::outerInternal:
                {
                    const labelHashSet cellSet(zoneCells);

                    const label nInt = mesh_.nInternalFaces();
                    const labelList& own = mesh_.faceOwner();
                    const labelList& nei = mesh_.faceNeighbour();
                    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

                    // Check all internal faces
                    for (label facei = 0; facei < nInt; facei++)
                    {
                        const bool inOwn = cellSet.found(own[facei]);
                        const bool inNei = cellSet.found(nei[facei]);

                        if (inOwn && !inNei)
                        {
                            selectedFaces[facei] = true;
                            flipMap[facei] = false;
                        }
                        else if (!inOwn && inNei)
                        {
                            selectedFaces[facei] = true;
                            flipMap[facei] = true;
                        }
                    }


                    // Get coupled cell status
                    boolList neiInSet(mesh_.nFaces()-nInt, false);

                    forAll(patches, patchi)
                    {
                        const polyPatch& pp = patches[patchi];

                        if (pp.coupled())
                        {
                            label facei = pp.start();
                            forAll(pp, i)
                            {
                                neiInSet[facei-nInt] =
                                    cellSet.found(own[facei]);
                                facei++;
                            }
                        }
                        else if
                        (
                            cellFaces_ != cellFaces::outerInternal
                        )
                        {
                            label facei = pp.start();
                            forAll(pp, i)
                            {
                                if (cellSet.found(own[facei]))
                                {
                                    selectedFaces[facei] = true;
                                }

                                facei++;
                            }
                        }
                    }

                    syncTools::swapBoundaryFaceList(mesh_, neiInSet);

                    // Check all boundary faces
                    forAll(patches, patchi)
                    {
                        const polyPatch& pp = patches[patchi];

                        if (pp.coupled())
                        {
                            label facei = pp.start();
                            forAll(pp, i)
                            {
                                const bool inOwn = cellSet.found(own[facei]);
                                const bool inNei = neiInSet[facei-nInt];

                                if (inOwn && !inNei)
                                {
                                    selectedFaces[facei] = true;
                                    flipMap[facei] = false;
                                }
                                else if (!inOwn && inNei)
                                {
                                    selectedFaces[facei] = true;
                                    flipMap[facei] = true;
                                }

                                facei++;
                            }
                        }
                    }
                    break;
                }
            }
        }

        if (zs.fZone.valid() && zs.fZone().name() != zoneName_)
        {
            const labelList& zoneFaces = zs.fZone();

            forAll(zoneFaces, zfi)
            {
                const face& f = mesh_.faces()[zoneFaces[zfi]];

                forAll(f, fp)
                {
                    selectedFaces[f[fp]] = true;
                }
            }
        }
    }

    moveUpdate_ = zoneGenerators_.moveUpdate();

    const labelList faceIndices(indices(selectedFaces));

    return zoneSet
    (
        new faceZone
        (
            zoneName_,
            faceIndices,
            boolList(flipMap, faceIndices),
            mesh_.faceZones(),
            moveUpdate_,
            true
        )
    );
}


// ************************************************************************* //
