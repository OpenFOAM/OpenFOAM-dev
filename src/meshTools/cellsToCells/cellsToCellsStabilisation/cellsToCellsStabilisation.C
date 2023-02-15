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

#include "cellsToCellsStabilisation.H"
#include "syncTools.H"
#include "wallPoint.H"
#include "WallLocationData.H"
#include "WallInfo.H"
#include "FaceCellWave.H"
#include "globalIndex.H"
#include "OBJstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellsToCellsStabilisation, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellsToCellsStabilisation::cellsToCellsStabilisation()
:
    stabilisation_(false),
    localStabilisationCells_(),
    stabilisationMapPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellsToCellsStabilisation::~cellsToCellsStabilisation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::cellsToCellsStabilisation::update
(
    const polyMesh& mesh,
    const PackedBoolList& cellCoupleds
)
{
    // Determine whether or not stabilisation is necessary
    stabilisation_ = false;
    forAll(cellCoupleds, celli)
    {
        if (!cellCoupleds[celli])
        {
            stabilisation_ = true;
            break;
        }
    }
    reduce(stabilisation_, orOp<bool>());

    // Quick return if nothing is to be done
    if (!stabilisation_) return;

    // Get some information regarding the cells on the other side of couplings
    List<remote> bFaceNbrProcCells(mesh.nFaces() - mesh.nInternalFaces());
    boolList bFaceNbrIsCoupled(mesh.nFaces() - mesh.nInternalFaces());
    forAll(bFaceNbrIsCoupled, bFacei)
    {
        const label owni = mesh.faceOwner()[bFacei + mesh.nInternalFaces()];
        bFaceNbrProcCells[bFacei] = remote(Pstream::myProcNo(), owni);
        bFaceNbrIsCoupled[bFacei] = cellCoupleds[owni];
    }
    syncTools::swapBoundaryFaceList(mesh, bFaceNbrProcCells);
    syncTools::swapBoundaryFaceList(mesh, bFaceNbrIsCoupled);

    // Determine the "cut" faces that separate coupled and un-coupled cellS
    typedef WallInfo<WallLocationData<wallPoint, remote>> info;
    DynamicList<label> cutFaces;
    DynamicList<info> cutFaceInfos;
    for (label facei = 0; facei < mesh.nFaces(); ++ facei)
    {
        const label owni = mesh.faceOwner()[facei];
        const bool ownIsCoupled = cellCoupleds[owni];

        if (facei < mesh.nInternalFaces())
        {
            const label nbri = mesh.faceNeighbour()[facei];
            const bool nbrIsCoupled = cellCoupleds[nbri];

            if (ownIsCoupled != nbrIsCoupled)
            {
                const label celli = ownIsCoupled ? owni : nbri;

                cutFaces.append(facei);
                cutFaceInfos.append
                (
                    info
                    (
                        remote(Pstream::myProcNo(), celli),
                        mesh.faceCentres()[facei],
                        0
                    )
                );
            }
        }
        else
        {
            const label bFacei = facei - mesh.nInternalFaces();

            const bool nbrIsCoupled = bFaceNbrIsCoupled[bFacei];

            if (!ownIsCoupled && nbrIsCoupled)
            {
                cutFaces.append(facei);
                cutFaceInfos.append
                (
                    info
                    (
                        bFaceNbrProcCells[bFacei],
                        mesh.faceCentres()[facei],
                        0
                    )
                );
            }
        }
    }

    // Wave the information about the cut faces' connected coupled cells into
    // the un-coupled cells. Base this wave on distance to the cut face.
    // Initialise coupled cells to have a distance of zero, so that we do not
    // waste time waving into coupled regions of the mesh.
    List<info> faceInfos(mesh.nFaces()), cellInfos(mesh.nCells());
    forAll(cellCoupleds, celli)
    {
        if (cellCoupleds[celli])
        {
            cellInfos[celli] =
                info
                (
                    remote(Pstream::myProcNo(), celli),
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
    forAll(cellCoupleds, celli)
    {
        if (!cellCoupleds[celli] && !cellInfos[celli].valid(wave.data()))
        {
            FatalErrorInFunction
                << "Un-coupled cell " << celli << " of mesh " << mesh.name()
                << " on processor " << Pstream::myProcNo() << " with centre "
                << "at " << mesh.cellCentres()[celli] << " was not connected "
                << "to a coupled cell by the stabilisation wave. This "
                << "indicates that an entire non-contiguous region of mesh "
                << "lies outside of the other mesh being mapped to. This is "
                << "not recoverable." << exit(FatalError);
        }
    }

    // Construct the cell to local stabilisation cell map
    const globalIndex cellGlobalIndex(mesh.nCells());
    localStabilisationCells_.resize(mesh.nCells());
    forAll(cellCoupleds, celli)
    {
        const remote& r = cellInfos[celli].data();
        localStabilisationCells_[celli] =
            cellCoupleds[celli]
          ? cellGlobalIndex.toGlobal(celli)
          : cellGlobalIndex.toGlobal(r.proci, r.elementi);
    }

    // Construct the distribution map, if necessary
    if (Pstream::parRun())
    {
        List<Map<label>> compactMap;
        stabilisationMapPtr_.reset
        (
            new distributionMap
            (
                cellGlobalIndex,
                localStabilisationCells_,
                compactMap
            )
        );
    }

    // Write out stabilisation connections
    if (debug)
    {
        OBJstream obj
        (
            typeName + "_" + mesh.name()
          + (Pstream::parRun() ? "_proc" + name(Pstream::myProcNo()) : "")
          + "_connections.obj"
        );

        pointField ccs(mesh.cellCentres());
        stabilise(ccs);

        forAll(ccs, celli)
        {
            const point& c = mesh.cellCentres()[celli];
            if (magSqr(c - ccs[celli]) == 0) continue;
            obj.write(linePointRef(ccs[celli], c));
        }
    }
}


// ************************************************************************* //
