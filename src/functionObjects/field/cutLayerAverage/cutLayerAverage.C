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

#include "cutLayerAverage.H"
#include "cutPolyIntegral.H"
#include "OSspecific.H"
#include "volPointInterpolation.H"
#include "writeFile.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cutLayerAverage, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        cutLayerAverage,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

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

    const scalarField& plotXs_;

    PtrList<DynamicList<label>> faceis_;

    PtrList<boolList> haveFaces_;

    inline label activeLayeri(const label layeri) const
    {
        return layeri % faceis_.size();
    }

public:

    faceData
    (
        const fvMesh& mesh,
        const label nActiveLayers,
        const scalarField& pointXs,
        const scalarField& plotXs
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
        plotXs_(plotXs),
        faceis_(nActiveLayers),
        haveFaces_(nActiveLayers)
    {
        forAll(faceis_, i)
        {
            faceis_.set(i, new DynamicList<label>(mesh.nFaces()));
            haveFaces_.set(i, new boolList(mesh.nFaces(), false));
        }
    }

    void clear(const label layeri)
    {
        UIndirectList<bool>
        (
            haveFaces_[activeLayeri(layeri)],
            faceis_[activeLayeri(layeri)]
        ) = false;

        faceis_[activeLayeri(layeri)].clear();
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
        const fvMesh& mesh,
        const label nActiveLayers,
        const scalarField& pointXs,
        const scalarField& plotXs
    )
    :
        faceData(mesh, nActiveLayers, pointXs, plotXs),
        faceCuts_(nActiveLayers),
        faceCutAreas_(nActiveLayers)
    {
        forAll(faceis_, i)
        {
            faceCuts_.set(i, new List<List<labelPair>>(mesh.nFaces()));
            faceCutAreas_.set(i, fieldPairPtr<vector>(mesh.nFaces()));
        }
    }

    void cache(const label celli, const label layeri)
    {
        forAll(cells_[celli], cellFacei)
        {
            const label facei = cells_[celli][cellFacei];

            if (haveFaces_[activeLayeri(layeri)][facei]) continue;

            faceis_[activeLayeri(layeri)].append(facei);

            haveFaces_[activeLayeri(layeri)][facei] = true;

            faceCuts_[activeLayeri(layeri)][facei] =
                cutPoly::faceCuts
                (
                    faces_[facei],
                    pointXs_,
                    plotXs_[layeri]
                );

            for (label sidei = 0; sidei <= 1; ++ sidei)
            {
                faceCutAreas_[activeLayeri(layeri)][sidei][facei] =
                    cutPoly::faceCutArea
                    (
                        faces_[facei],
                        faceAreas_[facei],
                        faceCuts_[activeLayeri(layeri)][facei],
                        points_,
                        pointXs_,
                        plotXs_[layeri],
                        sidei == 0
                    );
            }
        }
    }

    const List<List<labelPair>>& faceCuts(const label layeri) const
    {
        return faceCuts_[activeLayeri(layeri)];
    }

    const Pair<vectorField>& faceCutAreas(const label layeri) const
    {
        return faceCutAreas_[activeLayeri(layeri)];
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
        const fvMesh& mesh,
        const label nActiveLayers,
        const scalarField& pointXs,
        const scalarField& plotXs
    )
    :
        faceData(mesh, nActiveLayers, pointXs, plotXs),
        fcd_(fcd),
        pointFs_(nActiveLayers),
        faceFs_(nActiveLayers),
        faceCutFs_(nActiveLayers)
    {
        forAll(faceis_, i)
        {
            pointFs_.set(i, fieldPairPtr<scalar>(mesh.nPoints()));
            faceFs_.set(i, fieldPairPtr<scalar>(mesh.nFaces()));
            faceCutFs_.set(i, fieldPairPairPtr<scalar>(mesh.nFaces()));
        }
    }

    void cache(const label celli, const label layeri)
    {
        forAll(cells_[celli], cellFacei)
        {
            const label facei = cells_[celli][cellFacei];

            if (haveFaces_[activeLayeri(layeri)][facei]) continue;

            faceis_[activeLayeri(layeri)].append(facei);

            haveFaces_[activeLayeri(layeri)][facei] = true;

            forAll(faces_[facei], facePointi)
            {
                const label pointi = faces_[facei][facePointi];

                if (layeri > 0)
                {
                    pointFs_[activeLayeri(layeri)][0][pointi] =
                        (pointXs_[pointi] - plotXs_[layeri - 1])
                       /(plotXs_[layeri] - plotXs_[layeri - 1]);
                }

                if (layeri < plotXs_.size() - 1)
                {
                    pointFs_[activeLayeri(layeri)][1][pointi] =
                        (plotXs_[layeri + 1] - pointXs_[pointi])
                       /(plotXs_[layeri + 1] - plotXs_[layeri]);
                }
            }

            for
            (
                label fi = (layeri > 0 ? 0 : 1);
                fi <= (layeri < plotXs_.size() - 1 ? 1 : 0);
                fi ++
            )
            {
                faceFs_[activeLayeri(layeri)][fi][facei] =
                    cutPoly::faceAreaAverage
                    (
                        faces_[facei],
                        points_,
                        pointFs_[activeLayeri(layeri)][fi]
                    ).second();
            }

            for
            (
                label fi = (layeri > 0 ? 0 : 1);
                fi <= (layeri < plotXs_.size() - 1 ? 1 : 0);
                fi ++
            )
            {
                for (label sidei = 0; sidei <= 1; ++ sidei)
                {
                    const label layerj =
                        fi == 0 && sidei == 0 ? layeri - 1
                      : fi == 1 && sidei == 1 ? layeri + 1
                      : layeri;

                    if (layerj < 0 || layerj > plotXs_.size() - 1) continue;

                    Tuple2<vector, vector> integral =
                        cutPoly::faceCutAreaIntegral
                        (
                            faces_[facei],
                            faceAreas_[facei],
                            faceFs_[activeLayeri(layeri)][fi][facei],
                            fcd_.faceCuts(layerj)[facei],
                            points_,
                            pointFs_[activeLayeri(layeri)][fi],
                            pointXs_,
                            plotXs_[layerj],
                            sidei == 0
                        );

                    const scalar magSqrArea = magSqr(integral.first());

                    faceCutFs_[activeLayeri(layeri)][fi][sidei][facei] =
                        magSqrArea > vSmall
                      ? integral.second() & integral.first()/magSqrArea
                      : faceFs_[activeLayeri(layeri)][fi][facei];
                }
            }
        }
    }

    const Pair<scalarField>& pointFs(const label layeri) const
    {
        return pointFs_[activeLayeri(layeri)];
    }

    const Pair<scalarField>& faceFs(const label layeri) const
    {
        return faceFs_[activeLayeri(layeri)];
    }

    const Pair<Pair<scalarField>>& faceCutFs(const label layeri) const
    {
        return faceCutFs_[activeLayeri(layeri)];
    }
};

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::functionObjects::cutLayerAverage::outputPath() const
{
    return
        time_.globalPath()
       /writeFile::outputPrefix
       /(mesh_.name() != polyMesh::defaultRegion ? mesh_.name() : word())
       /name()
       /time_.name();
}


Foam::List<Foam::functionObjects::cutLayerAverage::weight>
Foam::functionObjects::cutLayerAverage::calcNonInterpolatingWeights
(
    const scalarField& pointXs,
    const scalarField& cellMinXs,
    const scalarField& cellMaxXs,
    const labelList& cellMinOrder,
    const scalarField& plotXs,
    const bool normalise
) const
{
    const cellList& cells = mesh_.cells();
    const cellEdgeAddressingList& cAddrs = cellEdgeAddressingList::New(mesh_);
    const scalarField& cellVolumes = mesh_.cellVolumes();
    const faceList& faces = mesh_.faces();
    const vectorField& faceAreas = mesh_.faceAreas();
    const pointField& faceCentres = mesh_.faceCentres();
    const pointField& points = mesh_.points();

    // Determine the largest number of layers spanned by a single cell, and
    // therefore the number for which data must be simultaneously maintained
    label nActiveLayers = 0;
    {
        label layeri = 0;
        forAll(cellMinOrder, cellMinOrderi)
        {
            const label celli = cellMinOrder[cellMinOrderi];

            // Find the next relevant layer
            while (cellMinXs[celli] > plotXs[layeri + 1]) layeri ++;

            // Find the first irrelevant layer
            label layerj = layeri;
            while (cellMaxXs[celli] > plotXs[layerj]) layerj ++;

            nActiveLayers = max(nActiveLayers, layerj - layeri);
        }
    }

    // Storage for face cut data
    faceCutData fcd(mesh_, nActiveLayers, pointXs, plotXs);

    // Generate weights for each cell in turn
    DynamicList<weight> dynWeights(cells.size()*2);
    label layeri = 0;
    forAll(cellMinOrder, cellMinOrderi)
    {
        const label celli = cellMinOrder[cellMinOrderi];

        // Find the next relevant layer and remove all data relating to
        // layers now behind the spans of the remaining cells
        while (cellMinXs[celli] > plotXs[layeri + 1])
        {
            fcd.clear(layeri);
            layeri ++;
        }

        // Loop over all relevant layer intervals
        label layerj = layeri;
        while (cellMaxXs[celli] > plotXs[layerj])
        {
            // Compute the face data as necessary
            fcd.cache(celli, layerj);
            fcd.cache(celli, layerj + 1);

            // Add a new weight
            dynWeights.append({celli, layerj, cellVolumes[celli]});

            // Left interval
            if (cellMinXs[celli] < plotXs[layerj])
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
                            fcd.faceCuts(layerj),
                            pointXs,
                            plotXs[layerj]
                        ),
                        faces,
                        faceAreas,
                        faceCentres,
                        fcd.faceCutAreas(layerj)[0],
                        points,
                        pointXs,
                        plotXs[layerj],
                        true
                    );
            }

            // Right interval
            if (cellMaxXs[celli] > plotXs[layerj + 1])
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
                            fcd.faceCuts(layerj + 1),
                            pointXs,
                            plotXs[layerj + 1]
                        ),
                        faces,
                        faceAreas,
                        faceCentres,
                        fcd.faceCutAreas(layerj + 1)[1],
                        points,
                        pointXs,
                        plotXs[layerj + 1],
                        false
                    );
            }

            layerj ++;
        }
    }

    // Transfer to non-dynamic storage
    List<weight> weights;
    weights.transfer(dynWeights);

    // Normalise, if requested
    if (normalise)
    {
        scalarField layerWeightSums(nLayers_, scalar(0));
        forAll(weights, weighti)
        {
            layerWeightSums[weights[weighti].layeri] += weights[weighti].value;
        }

        Pstream::listCombineGather(layerWeightSums, plusEqOp<scalar>());
        Pstream::listCombineScatter(layerWeightSums);

        forAll(weights, weighti)
        {
            weights[weighti].value /=
                (layerWeightSums[weights[weighti].layeri] + vSmall);
        }
    }

    return weights;
}


Foam::List<Foam::functionObjects::cutLayerAverage::weight>
Foam::functionObjects::cutLayerAverage::calcInterpolatingWeights
(
    const scalarField& pointXs,
    const scalarField& cellMinXs,
    const scalarField& cellMaxXs,
    const labelList& cellMinOrder,
    const scalarField& plotXs,
    const bool normalise
) const
{
    const cellList& cells = mesh_.cells();
    const cellEdgeAddressingList& cAddrs = cellEdgeAddressingList::New(mesh_);
    const scalarField& cellVolumes = mesh_.cellVolumes();
    const faceList& faces = mesh_.faces();
    const vectorField& faceAreas = mesh_.faceAreas();
    const pointField& faceCentres = mesh_.faceCentres();
    const pointField& points = mesh_.points();

    // Determine the largest number of layers spanned by a single cell, and
    // therefore the number for which data must be simultaneously maintained
    label nActiveLayers = 0;
    {
        label layeri = 0;
        forAll(cellMinOrder, cellMinOrderi)
        {
            const label celli = cellMinOrder[cellMinOrderi];

            // Find the next relevant layer
            while (cellMinXs[celli] > plotXs[layeri + 1]) layeri ++;

            // Find the first irrelevant layer
            label layerj = layeri;
            while (cellMaxXs[celli] > plotXs[max(layerj - 1, 0)]) layerj ++;

            nActiveLayers = max(nActiveLayers, layerj - layeri + 1);
        }
    }

    // Storage for face cut and face functions data
    faceCutData fcd(mesh_, nActiveLayers, pointXs, plotXs);
    faceFsData ffs(fcd, mesh_, nActiveLayers, pointXs, plotXs);

    // Generate weights for each cell in turn
    DynamicList<weight> dynWeights(cells.size()*2);
    label layeri = 0;
    forAll(cellMinOrder, cellMinOrderi)
    {
        const label celli = cellMinOrder[cellMinOrderi];

        // Find the next relevant layer and remove all data relating to
        // layers now behind the spans of the remaining cells
        while (cellMinXs[celli] > plotXs[layeri + 1])
        {
            if (layeri != 0) fcd.clear(layeri - 1);
            ffs.clear(layeri);
            layeri ++;
        }

        // Loop over all relevant layers
        label layerj = layeri;
        while (cellMaxXs[celli] > plotXs[max(layerj - 1, 0)])
        {
            // Compute the connected face data as necessary
            if (layerj != 0) fcd.cache(celli, layerj - 1);
            fcd.cache(celli, layerj);
            if (layerj != nLayers_ - 1) fcd.cache(celli, layerj + 1);
            ffs.cache(celli, layerj);

            // Add a new weight
            dynWeights.append({celli, layerj, 0});

            // Get the cell cuts for the middle of the layer. These will be
            // used twice. The cuts for the left and right will only be used
            // once and hence don't need to be stored. They are evaluated
            // inline below.
            const labelListList cCutsMid =
                cutPoly::cellCuts
                (
                    cells[celli],
                    cAddrs[celli],
                    faces,
                    fcd.faceCuts(layerj),
                    pointXs,
                    plotXs[layerj]
                );

            // Left interval
            if (layerj > 0 && cellMinXs[celli] < plotXs[layerj])
            {
                const scalar cellVF =
                    cutPoly::cellVolumeIntegral
                    (
                        cells[celli],
                        cAddrs[celli],
                        faceAreas,
                        faceCentres,
                        ffs.faceFs(layerj)[0]
                    ).second();

                // Add the whole cell's contribution
                dynWeights.last().value += cellVF;

                // Cut off anything before the left point
                if (cellMinXs[celli] < plotXs[layerj - 1])
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
                                fcd.faceCuts(layerj - 1),
                                pointXs,
                                plotXs[layerj - 1]
                            ),
                            faces,
                            faceAreas,
                            faceCentres,
                            ffs.faceFs(layerj)[0],
                            fcd.faceCutAreas(layerj - 1)[0],
                            ffs.faceCutFs(layerj)[0][0],
                            points,
                            ffs.pointFs(layerj)[0],
                            pointXs,
                            plotXs[layerj - 1],
                            true
                        ).second();
                }

                // Cut off anything after the middle point
                if (cellMaxXs[celli] > plotXs[layerj])
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
                            ffs.faceFs(layerj)[0],
                            fcd.faceCutAreas(layerj)[1],
                            ffs.faceCutFs(layerj)[0][1],
                            points,
                            ffs.pointFs(layerj)[0],
                            pointXs,
                            plotXs[layerj],
                            false
                        ).second();
                }
            }

            // Right interval
            if (layerj < nLayers_ - 1 && cellMaxXs[celli] > plotXs[layerj])
            {
                const scalar cellVF =
                    cutPoly::cellVolumeIntegral
                    (
                        cells[celli],
                        cAddrs[celli],
                        faceAreas,
                        faceCentres,
                        ffs.faceFs(layerj)[1]
                    ).second();

                // Add the whole cell's contribution
                dynWeights.last().value += cellVF;

                // Cut off anything before the middle point
                if (cellMinXs[celli] < plotXs[layerj])
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
                            ffs.faceFs(layerj)[1],
                            fcd.faceCutAreas(layerj)[0],
                            ffs.faceCutFs(layerj)[1][0],
                            points,
                            ffs.pointFs(layerj)[1],
                            pointXs,
                            plotXs[layerj],
                            true
                        ).second();
                }

                // Cut off anything after the right point
                if (cellMaxXs[celli] > plotXs[layerj + 1])
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
                                fcd.faceCuts(layerj + 1),
                                pointXs,
                                plotXs[layerj + 1]
                            ),
                            faces,
                            faceAreas,
                            faceCentres,
                            ffs.faceFs(layerj)[1],
                            fcd.faceCutAreas(layerj + 1)[1],
                            ffs.faceCutFs(layerj)[1][1],
                            points,
                            ffs.pointFs(layerj)[1],
                            pointXs,
                            plotXs[layerj + 1],
                            false
                        ).second();
                }
            }

            layerj ++;
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
        scalarField layerWeightSums(nLayers_, scalar(0));
        forAll(weights, weighti)
        {
            layerWeightSums[weights[weighti].layeri] += weights[weighti].value;
        }

        Pstream::listCombineGather(layerWeightSums, plusEqOp<scalar>());
        Pstream::listCombineScatter(layerWeightSums);

        forAll(weights, weighti)
        {
            weights[weighti].value /=
                (layerWeightSums[weights[weighti].layeri] + vSmall);
        }
    }
    else
    {
        forAll(weights, weighti)
        {
            if
            (
                weights[weighti].layeri == 0
             || weights[weighti].layeri == nLayers_ - 1
            )
            {
                weights[weighti].value *= 2;
            }
        }
    }

    return weights;
}


Foam::List<Foam::functionObjects::cutLayerAverage::weight>
Foam::functionObjects::cutLayerAverage::calcWeights
(
    const scalarField& pointXs,
    const scalarField& cellMinXs,
    const scalarField& cellMaxXs,
    const labelList& cellMinOrder,
    const scalarField& plotXs,
    const bool normalise
) const
{
    return
        interpolate_
      ? calcInterpolatingWeights
        (
            pointXs,
            cellMinXs,
            cellMaxXs,
            cellMinOrder,
            plotXs,
            normalise
        )
      : calcNonInterpolatingWeights
        (
            pointXs,
            cellMinXs,
            cellMaxXs,
            cellMinOrder,
            plotXs,
            normalise
        );
}


void Foam::functionObjects::cutLayerAverage::calcWeights()
{
    const cellList& cells = mesh_.cells();
    const faceList& faces = mesh_.faces();
    const pointField& points = mesh_.points();

    // If interpolating, then the layers and the plot points are coincident. If
    // not interpolating, then the layers lie in between the plot points, so
    // there is one more point than there are layers.
    const label nPlot = interpolate_ ? nLayers_ : nLayers_ + 1;

    // Calculate or get the point coordinates
    tmp<pointScalarField> tdistance =
        distanceName_ == word::null
      ? tmp<pointScalarField>(nullptr)
      : mesh_.foundObject<volScalarField>(distanceName_)
      ? volPointInterpolation::New(mesh_).interpolate
        (
            mesh_.lookupObject<volScalarField>(distanceName_)
        )
      : tmp<pointScalarField>
        (
            mesh_.lookupObject<pointScalarField>(distanceName_)
        );
    tmp<scalarField> tpointXs =
        distanceName_ == word::null
      ? points & direction_
      : tmp<scalarField>(tdistance().primitiveField());
    const scalarField& pointXs = tpointXs();

    // Determine face min and max coordinates
    scalarField cellMinXs(cells.size(), vGreat);
    scalarField cellMaxXs(cells.size(), -vGreat);
    forAll(cells, celli)
    {
        forAll(cells[celli], cellFacei)
        {
            const label facei = cells[celli][cellFacei];
            forAll(faces[facei], facePointi)
            {
                const label pointi = faces[facei][facePointi];
                cellMinXs[celli] = min(cellMinXs[celli], pointXs[pointi]);
                cellMaxXs[celli] = max(cellMaxXs[celli], pointXs[pointi]);
            }
        }
    }

    // Create orderings of the cells based on their min and max coordinates
    labelList cellMinOrder(cells.size());
    sortedOrder(cellMinXs, cellMinOrder);

    // Assume equal spacing to begin with
    const scalar xMin = gMin(pointXs), xMax = gMax(pointXs);
    scalarField plotXs
    (
        (xMin + scalarList(identityMap(nPlot))/(nPlot - 1)*(xMax - xMin))
    );

    // Names and fields for debug output of the counts, to observe the effect
    // of iterative improvement of the spacing
    wordList fieldNames;
    #define DeclareTypeFieldValues(Type, nullArg) \
        PtrList<Field<Type>> Type##FieldValues;
    FOR_ALL_FIELD_TYPES(DeclareTypeFieldValues);
    #undef DeclareTypeFieldValues

    // Iteratively optimise the spacing between the plot points to achieve an
    // approximately equal number of data points in each interval
    for (label iteri = 0; iteri < nOptimiseIter_ + debug; ++ iteri)
    {
        // Determine the count of faces that contribute to each layer
        const List<weight> weights =
            calcWeights
            (
                pointXs,
                cellMinXs,
                cellMaxXs,
                cellMinOrder,
                plotXs,
                false
            );
        const scalarField layerCounts
        (
            applyWeights<scalar>(weights, (1/mesh_.V())())
        );

        if (debug)
        {
            const label nFields0 = (2 + !interpolate_)*iteri;
            const label nFields = (2 + !interpolate_)*(iteri + 1);

            fieldNames.resize(nFields);
            #define ResizeTypeFieldValues(Type, nullArg) \
                Type##FieldValues.resize(nFields);
            FOR_ALL_FIELD_TYPES(ResizeTypeFieldValues);
            #undef ResizeTypeFieldValues

            if (!interpolate_)
            {
                const SubField<scalar> distance0s(plotXs, nLayers_);
                const SubField<scalar> distance1s(plotXs, nLayers_, 1);

                fieldNames[nFields0] = "distance-" + Foam::name(iteri);
                scalarFieldValues.set(nFields0, (distance0s + distance1s)/2);

                fieldNames[nFields0 + 1] = "thickness-" + Foam::name(iteri);
                scalarFieldValues.set(nFields0 + 1, distance1s - distance0s);
            }
            else
            {
                fieldNames[nFields0] = "distance-" + Foam::name(iteri);
                scalarFieldValues.set(nFields0, new scalarField(plotXs));
            }

            fieldNames[nFields - 1] = "count-" + Foam::name(iteri);
            scalarFieldValues.set(nFields - 1, new scalarField(layerCounts));

            if (iteri == nOptimiseIter_) break;
        }

        // Do a cumulative sum of the layer counts across all plot points
        scalarField plotSumCounts(nPlot, 0);
        for (label ploti = 0; ploti < nPlot - 1; ++ ploti)
        {
            plotSumCounts[ploti + 1] =
                plotSumCounts[ploti]
              + (
                    interpolate_
                  ? (layerCounts[ploti + 1] + layerCounts[ploti])/2
                  : layerCounts[ploti]
                );
        }

        // Compute the desired count in each interval
        const scalar plotDeltaCount = plotSumCounts.last()/(nPlot - 1);

        // Compute the new spacing between the points
        scalarField plot0Xs(plotXs);
        plotXs = -vGreat;
        plotXs.first() = xMin;
        label ploti = 1;
        for (label ploti0 = 0; ploti0 < nPlot - 1; ++ ploti0)
        {
            while
            (
                ploti < nPlot
             && plotSumCounts[ploti0 + 1] > ploti*plotDeltaCount
            )
            {
                const scalar f =
                    (ploti*plotDeltaCount - plotSumCounts[ploti0])
                   /(plotSumCounts[ploti0 + 1] - plotSumCounts[ploti0]);

                plotXs[ploti] = (1 - f)*plot0Xs[ploti0] + f*plot0Xs[ploti0 + 1];

                ploti ++;
            }
        }
        plotXs.last() = xMax;
    }

    if (debug)
    {
        mkDir(outputPath());

        formatter_->write
        (
            outputPath(),
            typeName + "_count",
            coordSet(labelList(nLayers_, 1)),
            fieldNames
            #define TypeFieldValuesParameter(Type, nullArg) \
                , Type##FieldValues
            FOR_ALL_FIELD_TYPES(TypeFieldValuesParameter)
            #undef TypeFieldValuesParameter
        );
    }

    // Finally, calculate the actual normalised interpolation weights
    weights_.reset
    (
        new List<weight>
        (
            calcWeights
            (
                pointXs,
                cellMinXs,
                cellMaxXs,
                cellMinOrder,
                plotXs,
                true
            )
        )
    );

    // Calculate plot coordinates and widths
    if (interpolate_)
    {
        layerDistances_.reset(new scalarField(plotXs));
    }
    else
    {
        const SubField<scalar> distance0s(plotXs, nLayers_);
        const SubField<scalar> distance1s(plotXs, nLayers_, 1);
        layerDistances_.reset(((distance0s + distance1s)/2).ptr());
        layerThicknesses_.reset((distance1s - distance0s).ptr());
    }

    // Calculate the plot positions
    layerPositions_.reset
    (
        applyWeights<vector>(weights_, mesh_.C()).ptr()
    );

    if (debug)
    {
        const List<weight> weights =
            calcWeights
            (
                pointXs,
                cellMinXs,
                cellMaxXs,
                cellMinOrder,
                plotXs,
                false
            );

        volTensorField::Internal layers
        (
            IOobject
            (
                name() + ":layers",
                mesh_.time().name(),
                mesh()
            ),
            mesh(),
            dimensionedTensor(dimless, tensor::zero)
        );

        forAll(weights, weighti)
        {
            const weight& w = weights[weighti];

            layers[w.celli][w.layeri % tensor::nComponents] =
                w.value/mesh_.V()[w.celli];
        }

        Info<< name() << ": Writing " << layers.name() << endl;

        layers.write();
    }
}


template<class Type>
inline Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::cutLayerAverage::applyWeights
(
    const List<weight>& weights,
    const VolInternalField<Type>& cellValues
) const
{
    tmp<Field<Type>> tLayerValues(new Field<Type>(nLayers_, Zero));

    forAll(weights, weighti)
    {
        tLayerValues.ref()[weights[weighti].layeri] +=
            weights[weighti].value*cellValues[weights[weighti].celli];
    }

    Pstream::listCombineGather(tLayerValues.ref(), plusEqOp<Type>());
    Pstream::listCombineScatter(tLayerValues.ref());

    return tLayerValues;
}


void Foam::functionObjects::cutLayerAverage::clear()
{
    weights_.clear();
    layerDistances_.clear();
    layerThicknesses_.clear();
    layerPositions_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cutLayerAverage::cutLayerAverage
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cutLayerAverage::~cutLayerAverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::cutLayerAverage::read(const dictionary& dict)
{
    const bool haveDirection = dict.found("direction");
    const bool haveDistance = dict.found("distance");
    if (haveDirection == haveDistance)
    {
        FatalIOErrorInFunction(dict)
            << "keywords direction and distance both "
            << (haveDirection ? "" : "un") << "defined in "
            << "dictionary " << dict.name()
            << exit(FatalIOError);
    }
    else if (haveDirection)
    {
        direction_ = normalised(dict.lookup<vector>("direction"));
        distanceName_ = word::null;
    }
    else if (haveDistance)
    {
        direction_ = vector::nan;
        distanceName_ = dict.lookup<word>("distance");
    }

    nLayers_ = dict.lookup<label>("nPoints");

    interpolate_ = dict.lookupOrDefault<bool>("interpolate", false);

    fields_ = dict.lookup<wordList>("fields");

    axis_ =
        coordSet::axisTypeNames_
        [
            dict.lookupOrDefault<word>
            (
                "axis",
                coordSet::axisTypeNames_[coordSet::axisType::DEFAULT]
            )
        ];

    formatter_ = setWriter::New(dict.lookup("setFormat"), dict);

    nOptimiseIter_ = dict.lookupOrDefault("nOptimiseIter", 2);

    clear();

    return true;
}


Foam::wordList Foam::functionObjects::cutLayerAverage::fields() const
{
    wordList result(fields_);

    if (distanceName_ != word::null)
    {
        result.append(distanceName_);
    }

    return result;
}


bool Foam::functionObjects::cutLayerAverage::execute()
{
    return true;
}


bool Foam::functionObjects::cutLayerAverage::write()
{
    if (!weights_.valid())
    {
        calcWeights();
    }

    const bool writeThickness =
        !interpolate_
     && (
            axis_ == coordSet::axisType::DEFAULT
         || axis_ == coordSet::axisType::DISTANCE
        );

    // Create list of field names
    wordList fieldNames;
    if (writeThickness)
    {
        fieldNames.append("thickness");
    }
    forAll(fields_, fieldi)
    {
        if
        (
            false
            #define FoundTypeField(Type, nullArg) \
              || foundObject<VolField<Type>>(fields_[fieldi])
            FOR_ALL_FIELD_TYPES(FoundTypeField)
            #undef FoundTypeField
        )
        {
            fieldNames.append(fields_[fieldi]);
        }
        else
        {
            cannotFindObject(fields_[fieldi]);
        }
    }

    // Calculate the field values
    #define DeclareTypeFieldValues(Type, nullArg) \
        PtrList<Field<Type>> Type##FieldValues(fieldNames.size());
    FOR_ALL_FIELD_TYPES(DeclareTypeFieldValues);
    #undef DeclareTypeFieldValues
    if (writeThickness)
    {
        scalarFieldValues.set(0, new scalarField(layerThicknesses_));
    }
    for (label fieldi = writeThickness; fieldi < fieldNames.size(); ++ fieldi)
    {
        #define CollapseTypeFields(Type, GeoField)                             \
            if (mesh_.foundObject<GeoField<Type>>(fieldNames[fieldi]))        \
            {                                                                 \
                const GeoField<Type>& field =                                 \
                    mesh_.lookupObject<GeoField<Type>>(fieldNames[fieldi]);   \
                                                                              \
                Type##FieldValues.set                                         \
                (                                                             \
                    fieldi,                                                   \
                    applyWeights<Type>(weights_, field) \
                );                                                            \
            }
        FOR_ALL_FIELD_TYPES(CollapseTypeFields, VolField);
        FOR_ALL_FIELD_TYPES(CollapseTypeFields, VolInternalField);
        #undef CollapseTypeFields
    }

    // Write
    if (Pstream::master() && layerPositions_->size())
    {
        mkDir(outputPath());

        formatter_->write
        (
            outputPath(),
            typeName,
            coordSet
            (
                identityMap(layerPositions_->size()),
                word::null,
                layerPositions_,
                coordSet::axisTypeNames_[coordSet::axisType::DISTANCE],
                layerDistances_,
                coordSet::axisTypeNames_[axis_]
            ),
            fieldNames
            #define TypeFieldValuesParameter(Type, nullArg) , Type##FieldValues
            FOR_ALL_FIELD_TYPES(TypeFieldValuesParameter)
            #undef TypeFieldValuesParameter
        );
    }

    return true;
}


void Foam::functionObjects::cutLayerAverage::movePoints
(
    const polyMesh& mesh
)
{
    if (&mesh == &mesh_)
    {
        clear();
    }
}


void Foam::functionObjects::cutLayerAverage::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        clear();
    }
}


void Foam::functionObjects::cutLayerAverage::mapMesh
(
    const polyMeshMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        clear();
    }
}


void Foam::functionObjects::cutLayerAverage::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        clear();
    }
}


// ************************************************************************* //
