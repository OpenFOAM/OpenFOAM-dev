/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "patchCutLayerAverage.H"
#include "cutPolyIntegral.H"
#include "OSspecific.H"
#include "volFields.H"
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
    defineTypeNameAndDebug(patchCutLayerAverage, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        patchCutLayerAverage,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::List<Foam::functionObjects::patchCutLayerAverage::weight>
Foam::functionObjects::patchCutLayerAverage::calcNonInterpolatingWeights
(
    const scalarField& pointXs,
    const scalarField& faceMinXs,
    const scalarField& faceMaxXs,
    const labelList& faceMinOrder,
    const scalarField& plotXs,
    const bool normalise
) const
{
    const polyPatch& pp = mesh_.boundaryMesh()[patchName_];
    const faceList& faces = pp.localFaces();
    const vectorField::subField& faceAreas = pp.faceAreas();
    const vectorField& faceNormals = pp.faceNormals();
    const pointField& points = pp.localPoints();

    // One-value for every point. Needed for the cutPoly routines.
    const scalarField pointOnes(points.size(), scalar(1));

    // Generate weights for each face in turn
    DynamicList<weight> dynWeights(faces.size()*2);
    label layeri = 0;
    forAll(faceMinOrder, i)
    {
        const label facei = faceMinOrder[i];

        const scalar a = faceAreas[facei] & faceNormals[facei];

        // Find the next relevant layer
        while (faceMinXs[facei] > plotXs[layeri + 1])
        {
            layeri ++;
        }

        // Loop over all relevant layer intervals
        label layerj = layeri;
        while (faceMaxXs[facei] > plotXs[layerj])
        {
            const scalar x0 = plotXs[layerj], x1 = plotXs[layerj + 1];

            // Add a new weight
            dynWeights.append({facei, layerj, a});

            // Left interval
            if (faceMinXs[facei] < x0)
            {
                dynWeights.last().value -=
                    cutPoly::faceCutAreaIntegral
                    (
                        faces[facei],
                        faceAreas[facei],
                        scalar(1),
                        cutPoly::faceCuts(faces[facei], pointXs, x0),
                        points,
                        pointOnes,
                        pointXs,
                        x0,
                        true
                    ).first() & faceNormals[facei];
            }

            // Right interval
            if (faceMaxXs[facei] > x1)
            {
                dynWeights.last().value -=
                    cutPoly::faceCutAreaIntegral
                    (
                        faces[facei],
                        faceAreas[facei],
                        scalar(1),
                        cutPoly::faceCuts(faces[facei], pointXs, x1),
                        points,
                        pointOnes,
                        pointXs,
                        x1,
                        false
                    ).first() & faceNormals[facei];
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


Foam::List<Foam::functionObjects::patchCutLayerAverage::weight>
Foam::functionObjects::patchCutLayerAverage::calcInterpolatingWeights
(
    const scalarField& pointXs,
    const scalarField& faceMinXs,
    const scalarField& faceMaxXs,
    const labelList& faceMinOrder,
    const scalarField& plotXs,
    const bool normalise
) const
{
    const polyPatch& pp = mesh_.boundaryMesh()[patchName_];
    const faceList& faces = pp.localFaces();
    const vectorField::subField& faceAreas = pp.faceAreas();
    const vectorField& faceNormals = pp.faceNormals();
    const pointField& points = pp.localPoints();

    // Values of the (piecewise linear) interpolation basis function. Different
    // for every interval. Is constantly being set on a sub-set of the
    // currently relevant points in the loop below.
    scalarField pointFs(points.size(), scalar(0));

    // Generate weights for each face in turn
    DynamicList<weight> dynWeights(faces.size()*2);
    label layeri = 0;
    forAll(faceMinOrder, i)
    {
        const label facei = faceMinOrder[i];

        const scalar a = faceAreas[facei] & faceNormals[facei];

        // Find the next relevant plot point
        while (faceMinXs[facei] > plotXs[layeri + 1])
        {
            layeri ++;
        }

        // Loop over all relevant plot points
        label layerj = layeri;
        while (layerj == 0 || faceMaxXs[facei] > plotXs[layerj - 1])
        {
            const scalar xLeft = layerj > 0 ? plotXs[layerj - 1] : NaN;
            const scalar xMid = plotXs[layerj];
            const scalar xRight =
                layerj < nLayers_ - 1 ? plotXs[layerj + 1] : NaN;

            // Add a new weight
            dynWeights.append({facei, layerj, 0});

            // Left interval
            if (layerj > 0 && faceMinXs[facei] < xMid)
            {
                forAll(faces[facei], facePointi)
                {
                    const label pointi = faces[facei][facePointi];
                    pointFs[pointi] =
                        (pointXs[pointi] - xLeft)/(xMid - xLeft);
                }

                const scalar f = faces[facei].average(points, pointFs);

                // Add the whole face's contribution
                dynWeights.last().value += f*a;

                // Cut off anything before the left point
                if (faceMinXs[facei] < xLeft)
                {
                    dynWeights.last().value -=
                        cutPoly::faceCutAreaIntegral
                        (
                            faces[facei],
                            faceAreas[facei],
                            f,
                            cutPoly::faceCuts(faces[facei], pointXs, xLeft),
                            points,
                            pointFs,
                            pointXs,
                            xLeft,
                            true
                        ).second() & faceNormals[facei];
                }

                // Cut off anything after the middle point
                if (faceMaxXs[facei] > xMid)
                {
                    dynWeights.last().value -=
                        cutPoly::faceCutAreaIntegral
                        (
                            faces[facei],
                            faceAreas[facei],
                            f,
                            cutPoly::faceCuts(faces[facei], pointXs, xMid),
                            points,
                            pointFs,
                            pointXs,
                            xMid,
                            false
                        ).second() & faceNormals[facei];
                }
            }

            // Right interval
            if (layerj < nLayers_ - 1 && faceMaxXs[facei] > xMid)
            {
                forAll(faces[facei], facePointi)
                {
                    const label pointi = faces[facei][facePointi];
                    pointFs[pointi] =
                        (xRight - pointXs[pointi])/(xRight - xMid);
                }

                const scalar f = faces[facei].average(points, pointFs);

                // Add the whole face's contribution
                dynWeights.last().value += f*a;

                // Cut off anything before the middle point
                if (faceMinXs[facei] < xMid)
                {
                    dynWeights.last().value -=
                        cutPoly::faceCutAreaIntegral
                        (
                            faces[facei],
                            faceAreas[facei],
                            f,
                            cutPoly::faceCuts(faces[facei], pointXs, xMid),
                            points,
                            pointFs,
                            pointXs,
                            xMid,
                            true
                        ).second() & faceNormals[facei];
                }

                // Cut off anything after the right point
                if (faceMaxXs[facei] > xRight)
                {
                    dynWeights.last().value -=
                        cutPoly::faceCutAreaIntegral
                        (
                            faces[facei],
                            faceAreas[facei],
                            f,
                            cutPoly::faceCuts(faces[facei], pointXs, xRight),
                            points,
                            pointFs,
                            pointXs,
                            xRight,
                            false
                        ).second() & faceNormals[facei];
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


Foam::List<Foam::functionObjects::patchCutLayerAverage::weight>
Foam::functionObjects::patchCutLayerAverage::calcWeights
(
    const scalarField& pointXs,
    const scalarField& faceMinXs,
    const scalarField& faceMaxXs,
    const labelList& faceMinOrder,
    const scalarField& plotXs,
    const bool normalise
) const
{
    return
        interpolate_
      ? calcInterpolatingWeights
        (
            pointXs,
            faceMinXs,
            faceMaxXs,
            faceMinOrder,
            plotXs,
            normalise
        )
      : calcNonInterpolatingWeights
        (
            pointXs,
            faceMinXs,
            faceMaxXs,
            faceMinOrder,
            plotXs,
            normalise
        );
}


template<class Type>
inline Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::patchCutLayerAverage::applyWeights
(
    const List<weight>& weights,
    const Field<Type>& faceValues
) const
{
    tmp<Field<Type>> tLayerValues(new Field<Type>(nLayers_, Zero));

    forAll(weights, weighti)
    {
        tLayerValues.ref()[weights[weighti].layeri] +=
            weights[weighti].value*faceValues[weights[weighti].facei];
    }

    Pstream::listCombineGather(tLayerValues.ref(), plusEqOp<Type>());
    Pstream::listCombineScatter(tLayerValues.ref());

    return tLayerValues;
}


void Foam::functionObjects::patchCutLayerAverage::initialise()
{
    const polyPatch& pp = mesh_.boundaryMesh()[patchName_];
    const faceList& faces = pp.localFaces();
    const pointField& points = pp.localPoints();

    // If interpolating, then the layers and the plot points are coincident. If
    // not interpolating, then the layers lie in between the plot points, so
    // there is one more point than there are layers.
    const label nPlot = interpolate_ ? nLayers_ : nLayers_ + 1;

    // Calculate point coordinates
    const scalarField pointXs(points & direction_);

    // Determine face min and max coordinates
    scalarField faceMinXs(faces.size(), vGreat);
    scalarField faceMaxXs(faces.size(), -vGreat);
    forAll(faces, facei)
    {
        forAll(faces[facei], facePointi)
        {
            const label pointi = faces[facei][facePointi];
            faceMinXs[facei] = min(faceMinXs[facei], pointXs[pointi]);
            faceMaxXs[facei] = max(faceMaxXs[facei], pointXs[pointi]);
        }
    }

    // Create orderings of the faces based on their min and max coordinates
    labelList faceMinOrder(faces.size());
    sortedOrder(faceMinXs, faceMinOrder);

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
                faceMinXs,
                faceMaxXs,
                faceMinOrder,
                plotXs,
                false
            );
        const scalarField layerCounts
        (
            applyWeights(weights, (1/pp.magFaceAreas())())
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
            while (plotSumCounts[ploti0 + 1] > ploti*plotDeltaCount)
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
    weights_ =
        calcWeights
        (
            pointXs,
            faceMinXs,
            faceMaxXs,
            faceMinOrder,
            plotXs,
            true
        );

    // Calculate plot coordinates and widths
    if (interpolate_)
    {
        layerDistances_ = plotXs;
    }
    else
    {
        const SubField<scalar> distance0s(plotXs, nLayers_);
        const SubField<scalar> distance1s(plotXs, nLayers_, 1);
        layerDistances_ = (distance0s + distance1s)/2;
        layerThicknesses_.reset((distance1s - distance0s).ptr());
    }

    // Calculate the plot positions
    layerPositions_ = applyWeights(weights_, pointField(pp.faceCentres()));
}


Foam::fileName Foam::functionObjects::patchCutLayerAverage::outputPath() const
{
    return
        time_.globalPath()
       /writeFile::outputPrefix
       /(mesh_.name() != polyMesh::defaultRegion ? mesh_.name() : word())
       /name()
       /time_.name();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::patchCutLayerAverage::patchCutLayerAverage
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

Foam::functionObjects::patchCutLayerAverage::~patchCutLayerAverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::patchCutLayerAverage::read(const dictionary& dict)
{
    patchName_ = dict.lookup<word>("patch");
    direction_ = normalised(dict.lookup<vector>("direction"));
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

    initialise();

    return true;
}


Foam::wordList Foam::functionObjects::patchCutLayerAverage::fields() const
{
    return fields_;
}


bool Foam::functionObjects::patchCutLayerAverage::execute()
{
    return true;
}


bool Foam::functionObjects::patchCutLayerAverage::write()
{
    const polyPatch& pp = mesh_.boundaryMesh()[patchName_];

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
        #define CollapseTypeFields(Type, nullArg)                             \
            if (mesh_.foundObject<VolField<Type>>(fieldNames[fieldi]))        \
            {                                                                 \
                const VolField<Type>& field =                                 \
                    mesh_.lookupObject<VolField<Type>>(fieldNames[fieldi]);   \
                                                                              \
                Type##FieldValues.set                                         \
                (                                                             \
                    fieldi,                                                   \
                    applyWeights(weights_, field.boundaryField()[pp.index()]) \
                );                                                            \
            }
        FOR_ALL_FIELD_TYPES(CollapseTypeFields);
        #undef CollapseTypeFields
    }

    // Write
    if (Pstream::master() && layerPositions_.size())
    {
        mkDir(outputPath());

        formatter_->write
        (
            outputPath(),
            typeName,
            coordSet
            (
                identityMap(layerPositions_.size()),
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


void Foam::functionObjects::patchCutLayerAverage::movePoints
(
    const polyMesh& mesh
)
{
    if (&mesh == &mesh_)
    {
        initialise();
    }
}


void Foam::functionObjects::patchCutLayerAverage::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        initialise();
    }
}


void Foam::functionObjects::patchCutLayerAverage::mapMesh
(
    const polyMeshMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        initialise();
    }
}


void Foam::functionObjects::patchCutLayerAverage::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        initialise();
    }
}


// ************************************************************************* //
