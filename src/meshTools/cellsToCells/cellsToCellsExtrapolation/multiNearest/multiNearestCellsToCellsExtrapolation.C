/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "multiNearestCellsToCellsExtrapolation.H"
#include "globalIndex.H"
#include "patchRegionSplit.H"
#include "wallPoint.H"
#include "WallLocationData.H"
#include "WallInfo.H"
#include "FaceCellWave.H"
#include "distributionMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellsToCellsExtrapolations
{
    defineTypeNameAndDebug(multiNearest, 0);
    addToRunTimeSelectionTable(cellsToCellsExtrapolation, multiNearest, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::cellsToCellsExtrapolations::multiNearest::extrapolateType
(
    Field<Type>& fld
) const
{
    if (!extrapolation_) return;

    // Communicate remote field values as necessary
    tmp<Field<Type>> srcFld;
    if (Pstream::parRun())
    {
        srcFld = fld.clone();
        extrapolationMapPtr_->distribute(srcFld.ref());
    }
    else
    {
        srcFld = tmp<Field<Type>>(fld);
    }

    // Set the values in the uncoupled cells
    forAll(uncoupledCells_, uncoupledCelli)
    {
        const label celli = uncoupledCells_[uncoupledCelli];

        fld[celli] = Zero;

        forAll(uncoupledCellLocalCells_[uncoupledCelli], i)
        {
            fld[celli] +=
                uncoupledCellWeights_[uncoupledCelli][i]
               *srcFld()[uncoupledCellLocalCells_[uncoupledCelli][i]];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellsToCellsExtrapolations::multiNearest::multiNearest()
:
    cellsToCellsExtrapolation(),
    uncoupledCellLocalCells_(),
    uncoupledCellWeights_(),
    extrapolationMapPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellsToCellsExtrapolations::multiNearest::~multiNearest()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::cellsToCellsExtrapolations::multiNearest::update
(
    const polyMesh& mesh,
    const PackedBoolList& cellCoupleds
)
{
    cellsToCellsExtrapolation::update(cellCoupleds);

    // Quick return if nothing is to be done
    if (!extrapolation_) return;

    // Global cell addressing
    const globalIndex globalCellIndex(mesh.nCells());

    // Get the list of faces which separate regions which do and do not need
    // extrapolation. For each face, also get the global index of the adjacent
    // cell which is coupled. Group them into continuous zones.
    labelListList zoneCutFaces, zoneCutFaceCoupledCells;
    {
        labelList cutFaces, cutFacesCoupledCells;
        getCutFaces
        (
            mesh,
            cellCoupleds,
            globalCellIndex,
            cutFaces,
            cutFacesCoupledCells
        );

        const patchRegionSplit cutFaceZones(mesh, cutFaces);

        labelList zoneNCutFaces(cutFaceZones.nRegions(), 0);
        forAll(cutFaces, cutFacei)
        {
            const label zonei = cutFaceZones[cutFacei];
            zoneNCutFaces[zonei] ++;
        }

        zoneCutFaces.resize(cutFaceZones.nRegions());
        zoneCutFaceCoupledCells.resize(cutFaceZones.nRegions());
        forAll(zoneNCutFaces, zonei)
        {
            zoneCutFaces[zonei].resize(zoneNCutFaces[zonei]);
            zoneCutFaceCoupledCells[zonei].resize(zoneNCutFaces[zonei]);
        }

        zoneNCutFaces = 0;
        forAll(cutFaces, cutFacei)
        {
            const label zonei = cutFaceZones[cutFacei];
            zoneCutFaces[zonei][zoneNCutFaces[zonei]] =
                cutFaces[cutFacei];
            zoneCutFaceCoupledCells[zonei][zoneNCutFaces[zonei]] =
                cutFacesCoupledCells[cutFacei];
            zoneNCutFaces[zonei] ++;
        }
    }

    // Initialise wave data. Give coupled cells a distance of zero, so that the
    // wave does not waste time propagating into coupled sections of the mesh.
    typedef WallInfo<WallLocationData<wallPoint, label>> info;
    List<info> faceInfos(mesh.nFaces()), cellInfos(mesh.nCells());
    forAll(cellCoupleds, celli)
    {
        if (cellCoupleds[celli])
        {
            cellInfos[celli] =
                info
                (
                    globalCellIndex.toGlobal(celli),
                    mesh.cellCentres()[celli],
                    0
                );
        }
    }

    //- Wave from each zone in turn. For each uncoupled cell that gets visited,
    //  note the global index of the nearest cell and store a weight equal to
    //  one-divided-by the distance.
    List<DynamicList<label>> uncoupledCellLocalCellsDyn(uncoupledCells_.size());
    List<DynamicList<scalar>> uncoupledCellWeightsDyn(uncoupledCells_.size());
    forAll(zoneCutFaces, zonei)
    {
        const labelList& cutFaces = zoneCutFaces[zonei];
        const labelList& cutFacesCoupledCells = zoneCutFaceCoupledCells[zonei];

        List<info> cutFaceInfos(cutFaces.size());
        forAll(cutFaces, cutFacei)
        {
            cutFaceInfos[cutFacei] =
                info
                (
                    cutFacesCoupledCells[cutFacei],
                    mesh.faceCentres()[cutFaces[cutFacei]],
                    0
                );
        }

        FaceCellWave<info> wave
        (
            mesh,
            cutFaces,
            cutFaceInfos,
            faceInfos,
            cellInfos,
            mesh.globalData().nTotalCells() + 1 // max iterations
        );

        forAll(uncoupledCells_, uncoupledCelli)
        {
            const label celli = uncoupledCells_[uncoupledCelli];

            if (cellInfos[celli].valid(wave.data()))
            {
                uncoupledCellLocalCellsDyn[uncoupledCelli].append
                (
                    cellInfos[celli].data()
                );
                uncoupledCellWeightsDyn[uncoupledCelli].append
                (
                    1/sqrt(cellInfos[celli].distSqr())
                );
            }
        }

        // Reset the wave data for the next zone
        forAll(uncoupledCells_, uncoupledCelli)
        {
            const label celli = uncoupledCells_[uncoupledCelli];

            cellInfos[celli] = info();

            forAll(mesh.cells()[celli], cellFacei)
            {
                const label facei = mesh.cells()[celli][cellFacei];

                faceInfos[facei] = info();
            }
        }
    }

    // Transfer to class storage
    uncoupledCellLocalCells_.resize(uncoupledCells_.size());
    uncoupledCellWeights_.resize(uncoupledCells_.size());
    forAll(uncoupledCells_, uncoupledCelli)
    {
        uncoupledCellLocalCells_[uncoupledCelli].transfer
        (
            uncoupledCellLocalCellsDyn[uncoupledCelli]
        );
        uncoupledCellWeights_[uncoupledCelli].transfer
        (
            uncoupledCellWeightsDyn[uncoupledCelli]
        );
    }

    // Check that the wave connected to all uncoupled cells
    forAll(uncoupledCells_, uncoupledCelli)
    {
        if (uncoupledCellLocalCells_[uncoupledCelli].empty())
        {
            const label celli = uncoupledCells_[uncoupledCelli];

            FatalErrorInFunction
                << "Un-coupled cell " << celli << " of mesh " << mesh.name()
                << " on processor " << Pstream::myProcNo() << " with centre "
                << "at " << mesh.cellCentres()[celli] << " was not connected "
                << "to a coupled cell by the extrapolation wave. This "
                << "indicates that an entire non-contiguous region of mesh "
                << "lies outside of the other mesh being mapped to. This is "
                << "not recoverable." << exit(FatalError);
        }
    }

    // Normalise the weights
    forAll(uncoupledCells_, uncoupledCelli)
    {
        scalar w = 0;
        forAll(uncoupledCellWeights_[uncoupledCelli], i)
        {
            w += uncoupledCellWeights_[uncoupledCelli][i];
        }
        forAll(uncoupledCellWeights_[uncoupledCelli], i)
        {
            uncoupledCellWeights_[uncoupledCelli][i] /= w;
        }
    }

    // Construct the distribution map, if necessary
    if (Pstream::parRun())
    {
        List<Map<label>> compactMap;
        extrapolationMapPtr_.reset
        (
            new distributionMap
            (
                globalCellIndex,
                uncoupledCellLocalCells_,
                compactMap
            )
        );
    }
}


#define implementExtrapolateType(Type, nullArg)                                \
    void Foam::cellsToCellsExtrapolations::multiNearest::extrapolate           \
    (                                                                          \
        Field<Type>& fld                                                       \
    ) const                                                                    \
    {                                                                          \
        extrapolateType(fld);                                                  \
    }
FOR_ALL_FIELD_TYPES(implementExtrapolateType);
#undef implementExtrapolateType


// ************************************************************************* //
