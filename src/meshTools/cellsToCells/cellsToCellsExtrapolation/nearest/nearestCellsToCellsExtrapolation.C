/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2026 OpenFOAM Foundation
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

#include "nearestCellsToCellsExtrapolation.H"
#include "globalIndex.H"
#include "wallPoint.H"
#include "WallLocationData.H"
#include "WallInfo.H"
#include "FaceCellWave.H"
#include "distributionMap.H"
#include "OBJstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellsToCellsExtrapolations
{
    defineTypeNameAndDebug(nearest, 0);
    addToRunTimeSelectionTable(cellsToCellsExtrapolation, nearest, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::cellsToCellsExtrapolations::nearest::extrapolateType
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

        fld[celli] = srcFld()[uncoupledCellLocalCells_[uncoupledCelli]];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellsToCellsExtrapolations::nearest::nearest()
:
    cellsToCellsExtrapolation(),
    uncoupledCellLocalCells_(),
    extrapolationMapPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellsToCellsExtrapolations::nearest::~nearest()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::cellsToCellsExtrapolations::nearest::update
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
    // cell which is coupled.
    labelList cutFaces, cutFacesCoupledCells;
    getCutFaces
    (
        mesh,
        cellCoupleds,
        globalCellIndex,
        cutFaces,
        cutFacesCoupledCells
    );

    // Wave the information about the cut faces' connected coupled cells into
    // the un-coupled cells. Base this wave on distance to the cut face.
    // Initialise coupled cells to have a distance of zero, so that we do not
    // waste time waving into coupled regions of the mesh.
    typedef WallInfo<WallLocationData<wallPoint, label>> info;
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
    FaceCellWave<info> wave
    (
        mesh,
        cutFaces,
        cutFaceInfos,
        faceInfos,
        cellInfos,
        mesh.globalData().nTotalCells() + 1 // max iterations
    );

    // Check that the wave connected to all un-coupled cells
    forAll(uncoupledCells_, uncoupledCelli)
    {
        const label celli = uncoupledCells_[uncoupledCelli];

        if (!cellInfos[celli].valid(wave.data()))
        {
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

    // Construct the cell to local extrapolation cell map
    uncoupledCellLocalCells_.resize(uncoupledCells_.size());
    forAll(uncoupledCells_, uncoupledCelli)
    {
        const label celli = uncoupledCells_[uncoupledCelli];

        uncoupledCellLocalCells_[uncoupledCelli] = cellInfos[celli].data();
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

    // Write out connections
    if (debug)
    {
        OBJstream obj
        (
            typeName + "_" + mesh.name()
          + (Pstream::parRun() ? "_proc" + name(Pstream::myProcNo()) : "")
          + "_connections.obj"
        );

        pointField ccs(mesh.cellCentres());
        extrapolate(ccs);

        forAll(ccs, celli)
        {
            const point& c = mesh.cellCentres()[celli];
            if (magSqr(c - ccs[celli]) == 0) continue;
            obj.write(linePointRef(ccs[celli], c));
        }
    }
}


#define implementExtrapolateType(Type, nullArg)                                \
    void Foam::cellsToCellsExtrapolations::nearest::extrapolate                \
    (                                                                          \
        Field<Type>& fld                                                       \
    ) const                                                                    \
    {                                                                          \
        extrapolateType(fld);                                                  \
    }
FOR_ALL_FIELD_TYPES(implementExtrapolateType);
#undef implementExtrapolateType


// ************************************************************************* //
