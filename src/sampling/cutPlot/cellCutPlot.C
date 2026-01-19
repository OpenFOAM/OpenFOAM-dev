/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2026 OpenFOAM Foundation
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

#include "cellCutPlot.H"
#include "cutPolyIntegral.H"
#include "generatedCellZone.H"
#include "OSspecific.H"
#include "setWriter.H"
#include "SubField.H"
#include "volFields.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cellCutPlot
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Allocate and return a pair of fields
template<class Type>
autoPtr<Pair<Field<Type>>> fieldPairPtr(const label size)
{
    return
        autoPtr<Pair<Field<Type>>>
        (
            new Pair<Field<Type>>(Field<Type>(size), Field<Type>(size))
        );
}


//- Allocate and return a pair of a pair of fields
template<class Type>
autoPtr<Pair<Pair<Field<Type>>>> fieldPairPairPtr(const label size)
{
    return
        autoPtr<Pair<Pair<Field<Type>>>>
        (
            new Pair<Pair<Field<Type>>>
            (
                Pair<Field<Type>>(Field<Type>(size), Field<Type>(size)),
                Pair<Field<Type>>(Field<Type>(size), Field<Type>(size))
            )
        );
}


//- Base class for classes which manage incomplete sets of face data
class faceData
{
protected:

    const cellList& cells_;

    const cellEdgeAddressingList& cAddrs_;

    const scalarField& cellVolumes_;

    const faceList& faces_;

    const vectorField& faceAreas_;

    const pointField& faceCentres_;

    const pointField& points_;

    const scalarField& pointXs_;

    const scalarField& cutXs_;

    PtrList<DynamicList<label>> faceis_;

    PtrList<boolList> haveFaces_;

    inline label activeCuti(const label cuti) const
    {
        return cuti % faceis_.size();
    }

public:

    faceData
    (
        const polyMesh& mesh,
        const label nActiveCuts,
        const scalarField& pointXs,
        const scalarField& cutXs
    )
    :
        cells_(mesh.cells()),
        cAddrs_(cellEdgeAddressingList::New(mesh)),
        cellVolumes_(mesh.cellVolumes()),
        faces_(mesh.faces()),
        faceAreas_(mesh.faceAreas()),
        faceCentres_(mesh.faceCentres()),
        points_(mesh.points()),
        pointXs_(pointXs),
        cutXs_(cutXs),
        faceis_(nActiveCuts),
        haveFaces_(nActiveCuts)
    {
        forAll(faceis_, i)
        {
            faceis_.set(i, new DynamicList<label>(mesh.nFaces()));
            haveFaces_.set(i, new boolList(mesh.nFaces(), false));
        }
    }

    void clear(const label cuti)
    {
        UIndirectList<bool>
        (
            haveFaces_[activeCuti(cuti)],
            faceis_[activeCuti(cuti)]
        ) = false;

        faceis_[activeCuti(cuti)].clear();
    }
};


//- Class to manage face cut data
class faceCutData
:
    public faceData
{
private:

    PtrList<List<List<labelPair>>> faceCuts_;

    PtrList<Pair<vectorField>> faceCutAreas_;

public:

    faceCutData
    (
        const polyMesh& mesh,
        const label nActiveCuts,
        const scalarField& pointXs,
        const scalarField& cutXs
    )
    :
        faceData(mesh, nActiveCuts, pointXs, cutXs),
        faceCuts_(nActiveCuts),
        faceCutAreas_(nActiveCuts)
    {
        forAll(faceis_, i)
        {
            faceCuts_.set(i, new List<List<labelPair>>(mesh.nFaces()));
            faceCutAreas_.set(i, fieldPairPtr<vector>(mesh.nFaces()));
        }
    }

    void cache(const label celli, const label cuti)
    {
        forAll(cells_[celli], cellFacei)
        {
            const label facei = cells_[celli][cellFacei];

            if (haveFaces_[activeCuti(cuti)][facei]) continue;

            faceis_[activeCuti(cuti)].append(facei);

            haveFaces_[activeCuti(cuti)][facei] = true;

            faceCuts_[activeCuti(cuti)][facei] =
                cutPoly::faceCuts
                (
                    faces_[facei],
                    pointXs_,
                    cutXs_[cuti]
                );

            for (label sidei = 0; sidei <= 1; ++ sidei)
            {
                faceCutAreas_[activeCuti(cuti)][sidei][facei] =
                    cutPoly::faceCutArea
                    (
                        faces_[facei],
                        faceAreas_[facei],
                        faceCuts_[activeCuti(cuti)][facei],
                        points_,
                        pointXs_,
                        cutXs_[cuti],
                        sidei == 0
                    );
            }
        }
    }

    const List<List<labelPair>>& faceCuts(const label cuti) const
    {
        return faceCuts_[activeCuti(cuti)];
    }

    const Pair<vectorField>& faceCutAreas(const label cuti) const
    {
        return faceCutAreas_[activeCuti(cuti)];
    }
};


//- Class to manage face basis function data
class faceFsData
:
    public faceData
{
private:

    const faceCutData& fcd_;

    PtrList<Pair<scalarField>> pointFs_;

    PtrList<Pair<scalarField>> faceFs_;

    PtrList<Pair<Pair<scalarField>>> faceCutFs_;


public:

    faceFsData
    (
        const faceCutData& fcd,
        const polyMesh& mesh,
        const label nActiveCuts,
        const scalarField& pointXs,
        const scalarField& cutXs
    )
    :
        faceData(mesh, nActiveCuts, pointXs, cutXs),
        fcd_(fcd),
        pointFs_(nActiveCuts),
        faceFs_(nActiveCuts),
        faceCutFs_(nActiveCuts)
    {
        forAll(faceis_, i)
        {
            pointFs_.set(i, fieldPairPtr<scalar>(mesh.nPoints()));
            faceFs_.set(i, fieldPairPtr<scalar>(mesh.nFaces()));
            faceCutFs_.set(i, fieldPairPairPtr<scalar>(mesh.nFaces()));
        }
    }

    void cache(const label celli, const label cuti)
    {
        forAll(cells_[celli], cellFacei)
        {
            const label facei = cells_[celli][cellFacei];

            if (haveFaces_[activeCuti(cuti)][facei]) continue;

            faceis_[activeCuti(cuti)].append(facei);

            haveFaces_[activeCuti(cuti)][facei] = true;

            forAll(faces_[facei], facePointi)
            {
                const label pointi = faces_[facei][facePointi];

                if (cuti > 0)
                {
                    pointFs_[activeCuti(cuti)][0][pointi] =
                        (pointXs_[pointi] - cutXs_[cuti - 1])
                       /(cutXs_[cuti] - cutXs_[cuti - 1]);
                }

                if (cuti < cutXs_.size() - 1)
                {
                    pointFs_[activeCuti(cuti)][1][pointi] =
                        (cutXs_[cuti + 1] - pointXs_[pointi])
                       /(cutXs_[cuti + 1] - cutXs_[cuti]);
                }
            }

            for
            (
                label fi = (cuti > 0 ? 0 : 1);
                fi <= (cuti < cutXs_.size() - 1 ? 1 : 0);
                fi ++
            )
            {
                faceFs_[activeCuti(cuti)][fi][facei] =
                    cutPoly::faceAreaAverage
                    (
                        faces_[facei],
                        points_,
                        pointFs_[activeCuti(cuti)][fi]
                    ).second();
            }

            for
            (
                label fi = (cuti > 0 ? 0 : 1);
                fi <= (cuti < cutXs_.size() - 1 ? 1 : 0);
                fi ++
            )
            {
                for (label sidei = 0; sidei <= 1; ++ sidei)
                {
                    const label cutj =
                        fi == 0 && sidei == 0 ? cuti - 1
                      : fi == 1 && sidei == 1 ? cuti + 1
                      : cuti;

                    if (cutj < 0 || cutj > cutXs_.size() - 1) continue;

                    Tuple2<vector, vector> integral =
                        cutPoly::faceCutAreaIntegral
                        (
                            faces_[facei],
                            faceAreas_[facei],
                            faceFs_[activeCuti(cuti)][fi][facei],
                            fcd_.faceCuts(cutj)[facei],
                            points_,
                            pointFs_[activeCuti(cuti)][fi],
                            pointXs_,
                            cutXs_[cutj],
                            sidei == 0
                        );

                    const scalar magSqrArea = magSqr(integral.first());

                    faceCutFs_[activeCuti(cuti)][fi][sidei][facei] =
                        magSqrArea > vSmall
                      ? integral.second() & integral.first()/magSqrArea
                      : faceFs_[activeCuti(cuti)][fi][facei];
                }
            }
        }
    }

    const Pair<scalarField>& pointFs(const label cuti) const
    {
        return pointFs_[activeCuti(cuti)];
    }

    const Pair<scalarField>& faceFs(const label cuti) const
    {
        return faceFs_[activeCuti(cuti)];
    }

    const Pair<Pair<scalarField>>& faceCutFs(const label cuti) const
    {
        return faceCutFs_[activeCuti(cuti)];
    }
};


List<weight> calcNonInterpolatingWeights
(
    const polyMesh& mesh,
    const generatedCellZone& zone,
    const scalarField& pointXs,
    const scalarField& zoneCellMinXs,
    const scalarField& zoneCellMaxXs,
    const labelList& zoneCellMinOrder,
    const scalarField& cutXs,
    const bool normalise
)
{
    const cellList& cells = mesh.cells();
    const cellEdgeAddressingList& cAddrs = cellEdgeAddressingList::New(mesh);
    const scalarField& cellVolumes = mesh.cellVolumes();
    const faceList& faces = mesh.faces();
    const vectorField& faceAreas = mesh.faceAreas();
    const pointField& faceCentres = mesh.faceCentres();
    const pointField& points = mesh.points();

    // Determine the largest number of cuts spanned by a single cell, and
    // therefore the number for which data must be simultaneously maintained
    label nActiveCuts = 0;
    {
        label cuti = 0;
        forAll(zoneCellMinOrder, zoneCellMinOrderi)
        {
            const label zoneCelli = zoneCellMinOrder[zoneCellMinOrderi];

            // Find the next relevant cut
            while (zoneCellMinXs[zoneCelli] > cutXs[cuti + 1]) cuti ++;

            // Find the first irrelevant cut
            label cutj = cuti;
            while (zoneCellMaxXs[zoneCelli] > cutXs[cutj]) cutj ++;

            nActiveCuts = max(nActiveCuts, cutj - cuti);
        }
    }

    // Storage for face cut data
    faceCutData fcd(mesh, nActiveCuts, pointXs, cutXs);

    // Generate weights for each cell in turn
    DynamicList<weight> dynWeights(zone.nCells()*2);
    label cuti = 0;
    forAll(zoneCellMinOrder, zoneCellMinOrderi)
    {
        const label zoneCelli = zoneCellMinOrder[zoneCellMinOrderi];
        const label celli = zone.celli(zoneCelli);

        // Find the next relevant cut and remove all data relating to
        // cuts now behind the spans of the remaining cells
        while
        (
            cuti < cutXs.size() - 1
         && zoneCellMinXs[zoneCelli] > cutXs[cuti + 1]
        )
        {
            fcd.clear(cuti);
            cuti ++;
        }

        // Loop over all relevant cuts
        label cutj = cuti;
        while
        (
            cutj < cutXs.size() - 1
         && zoneCellMaxXs[zoneCelli] > cutXs[cutj]
        )
        {
            // Compute the face data as necessary
            fcd.cache(celli, cutj);
            fcd.cache(celli, cutj + 1);

            // Add a new weight
            dynWeights.append({celli, cutj, cellVolumes[celli]});

            // Left interval
            if (zoneCellMinXs[zoneCelli] < cutXs[cutj])
            {
                dynWeights.last().value -=
                    cutPoly::cellCutVolume
                    (
                        cells[celli],
                        cAddrs[celli],
                        cellVolumes[celli],
                        cutPoly::cellCuts
                        (
                            cells[celli],
                            cAddrs[celli],
                            faces,
                            fcd.faceCuts(cutj),
                            pointXs,
                            cutXs[cutj]
                        ),
                        faces,
                        faceAreas,
                        faceCentres,
                        fcd.faceCutAreas(cutj)[0],
                        points,
                        pointXs,
                        cutXs[cutj],
                        true
                    );
            }

            // Right interval
            if (zoneCellMaxXs[zoneCelli] > cutXs[cutj + 1])
            {
                dynWeights.last().value -=
                    cutPoly::cellCutVolume
                    (
                        cells[celli],
                        cAddrs[celli],
                        cellVolumes[celli],
                        cutPoly::cellCuts
                        (
                            cells[celli],
                            cAddrs[celli],
                            faces,
                            fcd.faceCuts(cutj + 1),
                            pointXs,
                            cutXs[cutj + 1]
                        ),
                        faces,
                        faceAreas,
                        faceCentres,
                        fcd.faceCutAreas(cutj + 1)[1],
                        points,
                        pointXs,
                        cutXs[cutj + 1],
                        false
                    );
            }

            cutj ++;
        }
    }

    // Transfer to non-dynamic storage
    List<weight> weights;
    weights.transfer(dynWeights);

    // Normalise, if requested
    if (normalise)
    {
        scalarField cutWeightSums(cutXs.size() - 1, scalar(0));
        forAll(weights, weighti)
        {
            cutWeightSums[weights[weighti].cuti] += weights[weighti].value;
        }

        Pstream::listCombineGather(cutWeightSums, plusEqOp<scalar>());
        Pstream::listCombineScatter(cutWeightSums);

        forAll(weights, weighti)
        {
            weights[weighti].value /=
                (cutWeightSums[weights[weighti].cuti] + vSmall);
        }
    }

    return weights;
}


List<weight> calcInterpolatingWeights
(
    const polyMesh& mesh,
    const generatedCellZone& zone,
    const scalarField& pointXs,
    const scalarField& zoneCellMinXs,
    const scalarField& zoneCellMaxXs,
    const labelList& zoneCellMinOrder,
    const scalarField& cutXs,
    const bool normalise
)
{
    const cellList& cells = mesh.cells();
    const cellEdgeAddressingList& cAddrs = cellEdgeAddressingList::New(mesh);
    const scalarField& cellVolumes = mesh.cellVolumes();
    const faceList& faces = mesh.faces();
    const vectorField& faceAreas = mesh.faceAreas();
    const pointField& faceCentres = mesh.faceCentres();
    const pointField& points = mesh.points();

    // Determine the largest number of cuts spanned by a single cell, and
    // therefore the number for which data must be simultaneously maintained
    label nActiveCuts = 0;
    {
        label cuti = 0;
        forAll(zoneCellMinOrder, zoneCellMinOrderi)
        {
            const label zoneCelli = zoneCellMinOrder[zoneCellMinOrderi];

            // Find the next relevant cut
            while (zoneCellMinXs[zoneCelli] > cutXs[cuti + 1])
            {
                cuti ++;
            }

            // Find the first irrelevant cut
            label cutj = cuti;
            while (zoneCellMaxXs[zoneCelli] > cutXs[max(cutj - 1, 0)])
            {
                cutj ++;
            }

            nActiveCuts = max(nActiveCuts, cutj - cuti + 1);
        }
    }

    // Storage for face cut and face functions data
    faceCutData fcd(mesh, nActiveCuts, pointXs, cutXs);
    faceFsData ffs(fcd, mesh, nActiveCuts, pointXs, cutXs);

    // Generate weights for each cell in turn
    DynamicList<weight> dynWeights(zone.nCells()*2);
    label cuti = 0;
    forAll(zoneCellMinOrder, zoneCellMinOrderi)
    {
        const label zoneCelli = zoneCellMinOrder[zoneCellMinOrderi];
        const label celli = zone.celli(zoneCelli);

        // Find the next relevant cut and remove all data relating to
        // cuts now behind the spans of the remaining cells
        while
        (
            cuti < cutXs.size()
         && zoneCellMinXs[zoneCelli] > cutXs[cuti + 1]
        )
        {
            if (cuti != 0) fcd.clear(cuti - 1);
            ffs.clear(cuti);
            cuti ++;
        }

        // Loop over all relevant cuts
        label cutj = cuti;
        while
        (
            cutj < cutXs.size()
         && zoneCellMaxXs[zoneCelli] > cutXs[max(cutj - 1, 0)]
        )
        {
            // Compute the connected face data as necessary
            if (cutj != 0) fcd.cache(celli, cutj - 1);
            fcd.cache(celli, cutj);
            if (cutj != cutXs.size() - 1) fcd.cache(celli, cutj + 1);
            ffs.cache(celli, cutj);

            // Add a new weight
            dynWeights.append({celli, cutj, 0});

            // Get the cell cuts for the middle cut. These will be used twice.
            // The cuts for the left and right will only be used once and hence
            // don't need to be stored. They are evaluated inline below.
            const labelListList cCutsMid =
                cutPoly::cellCuts
                (
                    cells[celli],
                    cAddrs[celli],
                    faces,
                    fcd.faceCuts(cutj),
                    pointXs,
                    cutXs[cutj]
                );

            // Left interval
            if
            (
                cutj > 0
             && zoneCellMinXs[zoneCelli] < cutXs[cutj]
            )
            {
                const scalar cellVF =
                    cutPoly::cellVolumeIntegral
                    (
                        cells[celli],
                        cAddrs[celli],
                        faceAreas,
                        faceCentres,
                        ffs.faceFs(cutj)[0]
                    ).second();

                // Add the whole cell's contribution
                dynWeights.last().value += cellVF;

                // Cut off anything before the left point
                if (zoneCellMinXs[zoneCelli] < cutXs[cutj - 1])
                {
                    dynWeights.last().value -=
                        cutPoly::cellCutVolumeIntegral
                        (
                            cells[celli],
                            cAddrs[celli],
                            cellVolumes[celli],
                            cellVF/cellVolumes[celli],
                            cutPoly::cellCuts
                            (
                                cells[celli],
                                cAddrs[celli],
                                faces,
                                fcd.faceCuts(cutj - 1),
                                pointXs,
                                cutXs[cutj - 1]
                            ),
                            faces,
                            faceAreas,
                            faceCentres,
                            ffs.faceFs(cutj)[0],
                            fcd.faceCutAreas(cutj - 1)[0],
                            ffs.faceCutFs(cutj)[0][0],
                            points,
                            ffs.pointFs(cutj)[0],
                            pointXs,
                            cutXs[cutj - 1],
                            true
                        ).second();
                }

                // Cut off anything after the middle point
                if (zoneCellMaxXs[zoneCelli] > cutXs[cutj])
                {
                    dynWeights.last().value -=
                        cutPoly::cellCutVolumeIntegral
                        (
                            cells[celli],
                            cAddrs[celli],
                            cellVolumes[celli],
                            cellVF/cellVolumes[celli],
                            cCutsMid,
                            faces,
                            faceAreas,
                            faceCentres,
                            ffs.faceFs(cutj)[0],
                            fcd.faceCutAreas(cutj)[1],
                            ffs.faceCutFs(cutj)[0][1],
                            points,
                            ffs.pointFs(cutj)[0],
                            pointXs,
                            cutXs[cutj],
                            false
                        ).second();
                }
            }

            // Right interval
            if
            (
                cutj < cutXs.size() - 1
             && zoneCellMaxXs[zoneCelli] > cutXs[cutj]
            )
            {
                const scalar cellVF =
                    cutPoly::cellVolumeIntegral
                    (
                        cells[celli],
                        cAddrs[celli],
                        faceAreas,
                        faceCentres,
                        ffs.faceFs(cutj)[1]
                    ).second();

                // Add the whole cell's contribution
                dynWeights.last().value += cellVF;

                // Cut off anything before the middle point
                if (zoneCellMinXs[zoneCelli] < cutXs[cutj])
                {
                    dynWeights.last().value -=
                        cutPoly::cellCutVolumeIntegral
                        (
                            cells[celli],
                            cAddrs[celli],
                            cellVolumes[celli],
                            cellVF/cellVolumes[celli],
                            cCutsMid,
                            faces,
                            faceAreas,
                            faceCentres,
                            ffs.faceFs(cutj)[1],
                            fcd.faceCutAreas(cutj)[0],
                            ffs.faceCutFs(cutj)[1][0],
                            points,
                            ffs.pointFs(cutj)[1],
                            pointXs,
                            cutXs[cutj],
                            true
                        ).second();
                }

                // Cut off anything after the right point
                if (zoneCellMaxXs[zoneCelli] > cutXs[cutj + 1])
                {
                    dynWeights.last().value -=
                        cutPoly::cellCutVolumeIntegral
                        (
                            cells[celli],
                            cAddrs[celli],
                            cellVolumes[celli],
                            cellVF/cellVolumes[celli],
                            cutPoly::cellCuts
                            (
                                cells[celli],
                                cAddrs[celli],
                                faces,
                                fcd.faceCuts(cutj + 1),
                                pointXs,
                                cutXs[cutj + 1]
                            ),
                            faces,
                            faceAreas,
                            faceCentres,
                            ffs.faceFs(cutj)[1],
                            fcd.faceCutAreas(cutj + 1)[1],
                            ffs.faceCutFs(cutj)[1][1],
                            points,
                            ffs.pointFs(cutj)[1],
                            pointXs,
                            cutXs[cutj + 1],
                            false
                        ).second();
                }
            }

            cutj ++;
        }
    }

    // Transfer to non-dynamic storage
    List<weight> weights;
    weights.transfer(dynWeights);

    // Normalise, if requested. Otherwise, double the weight values on the ends
    // to account for the fact that these points only have half a basis function
    // contributing to their sums.
    if (normalise)
    {
        scalarField cutWeightSums(cutXs.size(), scalar(0));
        forAll(weights, weighti)
        {
            cutWeightSums[weights[weighti].cuti] += weights[weighti].value;
        }

        Pstream::listCombineGather(cutWeightSums, plusEqOp<scalar>());
        Pstream::listCombineScatter(cutWeightSums);

        forAll(weights, weighti)
        {
            weights[weighti].value /=
                (cutWeightSums[weights[weighti].cuti] + vSmall);
        }
    }
    else
    {
        forAll(weights, weighti)
        {
            if
            (
                weights[weighti].cuti == 0
             || weights[weighti].cuti == cutXs.size() - 1
            )
            {
                weights[weighti].value *= 2;
            }
        }
    }

    return weights;
}


List<weight> calcWeights
(
    const polyMesh& mesh,
    const generatedCellZone& zone,
    const scalarField& pointXs,
    const scalarField& cellMinXs,
    const scalarField& cellMaxXs,
    const labelList& cellMinOrder,
    const scalarField& cutXs,
    const bool interpolate,
    const bool normalise
)
{
    return
        (
            interpolate
          ? calcInterpolatingWeights
          : calcNonInterpolatingWeights
        )
        (
            mesh,
            zone,
            pointXs,
            cellMinXs,
            cellMaxXs,
            cellMinOrder,
            cutXs,
            normalise
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cellCutPlot
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::List<Foam::cellCutPlot::weight> Foam::cellCutPlot::calcWeights
(
    const polyMesh& mesh,
    const generatedCellZone& zone,
    const scalarField& pointXs,
    const scalarField& cutXs,
    const bool interpolate,
    const bool normalise
)
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();

    // Determine cell min and max coordinates
    scalarField zoneCellMinXs(zone.nCells(), vGreat);
    scalarField zoneCellMaxXs(zone.nCells(), -vGreat);
    forAll(zoneCellMinXs, zoneCelli)
    {
        const label celli = zone.celli(zoneCelli);
        forAll(cells[celli], cellFacei)
        {
            const label facei = cells[celli][cellFacei];
            forAll(faces[facei], facePointi)
            {
                const label pointi = faces[facei][facePointi];
                zoneCellMinXs[zoneCelli] =
                    min(zoneCellMinXs[zoneCelli], pointXs[pointi]);
                zoneCellMaxXs[zoneCelli] =
                    max(zoneCellMaxXs[zoneCelli], pointXs[pointi]);
            }
        }
    }

    // Create orderings of the cells based on their min and max coordinates
    labelList zoneCellMinOrder(zone.nCells());
    sortedOrder(zoneCellMinXs, zoneCellMinOrder);

    // Calculate and return the weights
    return
        calcWeights
        (
            mesh,
            zone,
            pointXs,
            zoneCellMinXs,
            zoneCellMaxXs,
            zoneCellMinOrder,
            cutXs,
            interpolate,
            normalise
        );
}


Foam::tmp<Foam::scalarField> Foam::cellCutPlot::calcCutXs
(
    const polyMesh& mesh,
    const generatedCellZone& zone,
    const scalarField& pointXs,
    const bool interpolate,
    const label nCuts,
    const label nIter,
    const bool debug,
    const word& functionName,
    const setWriter& functionFormatter
)
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();

    const label nIntervals = interpolate ? nCuts : nCuts - 1;

    // Determine face min and max coordinates
    scalarField zoneCellMinXs(zone.nCells(), vGreat);
    scalarField zoneCellMaxXs(zone.nCells(), -vGreat);
    forAll(zoneCellMinXs, zoneCelli)
    {
        const label celli = zone.celli(zoneCelli);
        forAll(cells[celli], cellFacei)
        {
            const label facei = cells[celli][cellFacei];
            forAll(faces[facei], facePointi)
            {
                const label pointi = faces[facei][facePointi];
                zoneCellMinXs[zoneCelli] =
                    min(zoneCellMinXs[zoneCelli], pointXs[pointi]);
                zoneCellMaxXs[zoneCelli] =
                    max(zoneCellMaxXs[zoneCelli], pointXs[pointi]);
            }
        }
    }

    // Create orderings of the cells based on their min and max coordinates
    labelList zoneCellMinOrder(zone.nCells());
    sortedOrder(zoneCellMinXs, zoneCellMinOrder);

    // Assume equal spacing to begin with
    scalar xMin = gMin(pointXs), xMax = gMax(pointXs);
    xMin -= max(rootVSmall, 2*small*mag(xMin));
    xMax += max(rootVSmall, 2*small*mag(xMax));
    tmp<scalarField> tcutXs =
        (xMin + scalarList(identityMap(nCuts))/(nCuts - 1)*(xMax - xMin));
    scalarField& cutXs = tcutXs.ref();
    cutXs.first() = xMin;
    cutXs.last() = xMax;

    // Names and fields for debug output of the counts, to observe the effect
    // of iterative improvement of the spacing
    wordList fieldNames;
    #define DeclareTypeFieldValues(Type, nullArg) \
        PtrList<Field<Type>> Type##FieldValues;
    FOR_ALL_FIELD_TYPES(DeclareTypeFieldValues);
    #undef DeclareTypeFieldValues

    // Iteratively optimise the spacing between the cut points to achieve an
    // approximately equal number of data points in each interval
    for (label iteri = 0; iteri < nIter + debug; ++ iteri)
    {
        // Determine the count of faces that contribute to each interval
        const List<weight> weights =
            calcWeights
            (
                mesh,
                zone,
                pointXs,
                zoneCellMinXs,
                zoneCellMaxXs,
                zoneCellMinOrder,
                cutXs,
                interpolate,
                false
            );
        const scalarField intervalCounts
        (
            cutPlot::applyWeights(nIntervals, weights, (1/mesh.cellVolumes())())
        );

        if (debug)
        {
            const label nFields0 = (2 + !interpolate)*iteri;
            const label nFields = (2 + !interpolate)*(iteri + 1);

            fieldNames.resize(nFields);
            #define ResizeTypeFieldValues(Type, nullArg) \
                Type##FieldValues.resize(nFields);
            FOR_ALL_FIELD_TYPES(ResizeTypeFieldValues);
            #undef ResizeTypeFieldValues

            if (!interpolate)
            {
                const SubField<scalar> distance0s(cutXs, nIntervals);
                const SubField<scalar> distance1s(cutXs, nIntervals, 1);

                fieldNames[nFields0] = "distance-" + Foam::name(iteri);
                scalarFieldValues.set(nFields0, (distance0s + distance1s)/2);

                fieldNames[nFields0 + 1] = "thickness-" + Foam::name(iteri);
                scalarFieldValues.set(nFields0 + 1, distance1s - distance0s);
            }
            else
            {
                fieldNames[nFields0] = "distance-" + Foam::name(iteri);
                scalarFieldValues.set(nFields0, new scalarField(cutXs));
            }

            fieldNames[nFields - 1] = "count-" + Foam::name(iteri);
            scalarFieldValues.set(nFields - 1, new scalarField(intervalCounts));

            if (iteri == nIter) break;
        }

        // Do a cumulative sum of the interval counts across all cut points
        scalarField cutSumCounts(nCuts, 0);
        for (label cuti = 0; cuti < nCuts - 1; ++ cuti)
        {
            cutSumCounts[cuti + 1] =
                cutSumCounts[cuti]
              + (
                    interpolate
                  ? (intervalCounts[cuti + 1] + intervalCounts[cuti])/2
                  : intervalCounts[cuti]
                );
        }

        // Compute the desired count in each interval
        const scalar intervalCount = cutSumCounts.last()/(nCuts - 1);

        // Compute the new spacing between the points
        scalarField cut0Xs(cutXs);
        cutXs = -vGreat;
        cutXs.first() = xMin;
        label cuti = 1;
        for (label cuti0 = 0; cuti0 < nCuts - 1; ++ cuti0)
        {
            while
            (
                cuti < nCuts
             && cutSumCounts[cuti0 + 1] > cuti*intervalCount
            )
            {
                const scalar f =
                    (cuti*intervalCount - cutSumCounts[cuti0])
                   /(cutSumCounts[cuti0 + 1] - cutSumCounts[cuti0]);

                cutXs[cuti] = (1 - f)*cut0Xs[cuti0] + f*cut0Xs[cuti0 + 1];

                cuti ++;
            }
        }
        cutXs.last() = xMax;
    }

    if (debug)
    {
        const fileName outputPath =
            mesh.time().globalPath()
           /functionObjects::writeFile::outputPrefix
           /(
               mesh.name() != polyMesh::defaultRegion
             ? mesh.name()
             : word()
            )
           /functionName
           /mesh.time().name();

        mkDir(outputPath);

        functionFormatter.write
        (
            outputPath,
            functionName + "_count",
            coordSet(labelList(nIntervals, 1)),
            fieldNames
            #define TypeFieldValuesParameter(Type, nullArg) \
                , Type##FieldValues
            FOR_ALL_FIELD_TYPES(TypeFieldValuesParameter)
            #undef TypeFieldValuesParameter
        );
    }

    return tcutXs;
}


void Foam::cellCutPlot::writeLayers
(
    const fvMesh& mesh,
    const List<cellCutPlot::weight>& weights,
    const word& functionName
)
{
    volTensorField::Internal layers
    (
        IOobject
        (
            functionName + ":layers",
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedTensor(dimless, tensor::uniform(-1))
    );

    forAll(weights, weighti)
    {
        const weight& w = weights[weighti];

        layers[w.elementi] = tensor::zero;
    }

    forAll(weights, weighti)
    {
        const weight& w = weights[weighti];

        layers[w.elementi][w.cuti % tensor::nComponents] =
            w.value/mesh.cellVolumes()[w.elementi];
    }

    Info<< functionName << ": Writing " << layers.name() << endl;

    layers.write();
}


// ************************************************************************* //
