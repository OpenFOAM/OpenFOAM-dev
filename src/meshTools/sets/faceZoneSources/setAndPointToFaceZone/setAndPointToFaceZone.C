/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "setAndPointToFaceZone.H"
#include "FaceCellWave.H"
#include "faceZoneSet.H"
#include "indexedOctree.H"
#include "minData.H"
#include "treeDataCell.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setAndPointToFaceZone, 0);
    addToRunTimeSelectionTable
    (
        topoSetSource,
        setAndPointToFaceZone,
        word
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setAndPointToFaceZone::setAndPointToFaceZone
(
    const polyMesh& mesh,
    const word& setName,
    const vector& point
)
:
    topoSetSource(mesh),
    setName_(setName),
    point_(point)
{}


Foam::setAndPointToFaceZone::setAndPointToFaceZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.lookup<word>("set")),
    point_(dict.lookup<vector>("point", dimLength))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::setAndPointToFaceZone::~setAndPointToFaceZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::setAndPointToFaceZone::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (!isA<faceZoneSet>(set))
    {
        WarningInFunction
            << "Operation only allowed on a faceZoneSet."
            << endl;
        return;
    }

    // Cast to get the zone
    faceZoneSet& fzSet = refCast<faceZoneSet>(set);

    // Load the set
    faceSet loadedSet(mesh_, setName_);

    // Allocate new topology
    DynamicList<label> newAddressing(fzSet.addressing().size());
    DynamicList<bool> newFlipMap(fzSet.flipMap().size());

    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all faces from faceSet " << setName_
            << " ..." << endl;

        // Construct a wave to transport a minimum index
        List<minData> allFaceData(mesh_.nFaces());
        List<minData> allCellData(mesh_.nCells());

        // Set the existing faces of the set to a higher index than that of
        // uninitialised data
        forAllConstIter(faceSet, loadedSet, iter)
        {
            allFaceData[iter.key()] = minData(-1);
        }

        // Find the cell containing the given point and give that an even
        // higher index
        const label celli = mesh_.cellTree().findInside(point_);
        if (returnReduce(celli, maxOp<label>()) == -1)
        {
            FatalErrorInFunction
                << "Point " << point_ << " was not found in the mesh"
                << exit(FatalError);
        }
        if (celli != -1)
        {
            allCellData[celli] = minData(0);
        }

        // Create seed faces as all those of the cell containing the given
        // point. Give these the same index as the cell.
        const label nSeedFaces = celli != -1 ? mesh_.cells()[celli].size() : 0;
        labelList seedFaces(nSeedFaces);
        List<minData> seedFaceData(nSeedFaces);
        if (celli != -1)
        {
            forAll(mesh_.cells()[celli], cellFacei)
            {
                const label facei = mesh_.cells()[celli][cellFacei];
                seedFaces[cellFacei] = facei;
                seedFaceData = minData(0);
            }
        }

        // Wave from the point to the set faces
        FaceCellWave<minData> wave
        (
            mesh_,
            seedFaces,
            seedFaceData,
            allFaceData,
            allCellData,
            mesh().globalData().nTotalCells()+1
        );

        // Initialise as a copy
        newAddressing = fzSet.addressing();
        newFlipMap = fzSet.flipMap();

        // Unpack. Each set face should be at the boundary between cells that
        // were and were not visited by the wave. If it is the neighbour that
        // was visited by the wave, then set the flip map.
        forAllConstIter(faceSet, loadedSet, iter)
        {
            const label facei = iter.key();

            const label owni = mesh_.faceOwner()[facei];
            const bool ownValid = allCellData[owni].valid(wave.data());

            if (mesh_.isInternalFace(facei))
            {
                const label nbri = mesh_.faceNeighbour()[facei];
                const bool nbrValid = allCellData[nbri].valid(wave.data());

                if (ownValid && nbrValid)
                {
                    FatalErrorInFunction
                        << "Internal face #" << facei << " at "
                        << mesh_.faceCentres()[facei]
                        << " was visited from both sides by a wave from "
                        << point_ << exit(FatalError);
                }

                if (!ownValid && !nbrValid)
                {
                    FatalErrorInFunction
                        << "Internal face #" << facei << " at "
                        << mesh_.faceCentres()[facei]
                        << " was not visited by a wave from "
                        << point_ << exit(FatalError);
                }
            }
            else
            {
                if (!ownValid)
                {
                    FatalErrorInFunction
                        << "Boundary face #" << facei << " at "
                        << mesh_.faceCentres()[facei]
                        << " was not visited by a wave from "
                        << point_ << exit(FatalError);
                }
            }

            newAddressing.append(facei);
            newFlipMap.append(!ownValid);
        }
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all faces from faceSet " << setName_
            << " ..." << endl;

        // Remove everything in the set from the zone. Don't have to worry
        // about computing the flips here, seeing as we are only removing
        // faces.
        DynamicList<label> newAddressing(fzSet.addressing().size());
        DynamicList<bool> newFlipMap(fzSet.flipMap().size());
        forAll(fzSet.addressing(), i)
        {
            if (!loadedSet.found(fzSet.addressing()[i]))
            {
                newAddressing.append(fzSet.addressing()[i]);
                newFlipMap.append(fzSet.flipMap()[i]);
            }
        }
    }

    // Transfer into the zone
    fzSet.addressing().transfer(newAddressing);
    fzSet.flipMap().transfer(newFlipMap);
    fzSet.updateSet();
}


// ************************************************************************* //
