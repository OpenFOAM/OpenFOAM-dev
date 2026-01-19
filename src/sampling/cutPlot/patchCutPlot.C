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

#include "patchCutPlot.H"
#include "cutPolyIntegral.H"
#include "OSspecific.H"
#include "setWriter.H"
#include "SubField.H"
#include "volFields.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace patchCutPlot
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

List<weight> calcNonInterpolatingWeights
(
    const faceList& faces,
    const UList<vector>& faceAreas,
    const UList<vector>& faceNormals,
    const pointField& points,
    const scalarField& pointXs,
    const scalarField& faceMinXs,
    const scalarField& faceMaxXs,
    const labelList& faceMinOrder,
    const scalarField& cutXs,
    const bool normalise
)
{
    // Generate weights for each face in turn
    DynamicList<weight> dynWeights(faces.size()*2);
    label cuti = 0;
    forAll(faceMinOrder, i)
    {
        const label facei = faceMinOrder[i];

        const scalar a = faceAreas[facei] & faceNormals[facei];

        // Find the next relevant cut
        while
        (
            cuti < cutXs.size() - 1
         && faceMinXs[facei] > cutXs[cuti + 1]
        )
        {
            cuti ++;
        }

        // Loop over all relevant cut intervals
        label cutj = cuti;
        while
        (
            cutj < cutXs.size() - 1
         && faceMaxXs[facei] > cutXs[cutj]
        )
        {
            // Add a new weight
            dynWeights.append({facei, cutj, a});

            // Left interval
            if (faceMinXs[facei] < cutXs[cutj])
            {
                dynWeights.last().value -=
                    cutPoly::faceCutArea
                    (
                        faces[facei],
                        faceAreas[facei],
                        cutPoly::faceCuts
                        (
                            faces[facei],
                            pointXs,
                            cutXs[cutj]
                        ),
                        points,
                        pointXs,
                        cutXs[cutj],
                        true
                    ) & faceNormals[facei];
            }

            // Right interval
            if (faceMaxXs[facei] > cutXs[cutj + 1])
            {
                dynWeights.last().value -=
                    cutPoly::faceCutArea
                    (
                        faces[facei],
                        faceAreas[facei],
                        cutPoly::faceCuts
                        (
                            faces[facei],
                            pointXs,
                            cutXs[cutj + 1]
                        ),
                        points,
                        pointXs,
                        cutXs[cutj + 1],
                        false
                    ) & faceNormals[facei];
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
    const faceList& faces,
    const UList<vector>& faceAreas,
    const UList<vector>& faceNormals,
    const pointField& points,
    const scalarField& pointXs,
    const scalarField& faceMinXs,
    const scalarField& faceMaxXs,
    const labelList& faceMinOrder,
    const scalarField& cutXs,
    const bool normalise
)
{
    // Values of the (piecewise linear) interpolation basis function. Different
    // for every interval. Is constantly being set on a sub-set of the
    // currently relevant points in the loop below.
    scalarField pointFs(points.size(), scalar(0));

    // Generate weights for each face in turn
    DynamicList<weight> dynWeights(faces.size()*2);
    label cuti = 0;
    forAll(faceMinOrder, i)
    {
        const label facei = faceMinOrder[i];

        const scalar a = faceAreas[facei] & faceNormals[facei];

        // Find the next relevant cut
        while
        (
            cuti < cutXs.size()
         && faceMinXs[facei] > cutXs[cuti + 1]
        )
        {
            cuti ++;
        }

        // Loop over all relevant cuts
        label cutj = cuti;
        while
        (
            cutj < cutXs.size()
         && faceMaxXs[facei] > cutXs[max(cutj - 1, 0)]
        )
        {
            // Add a new weight
            dynWeights.append({facei, cutj, 0});

            // Left interval
            if (cutj > 0 && faceMinXs[facei] < cutXs[cutj])
            {
                // Update the basis function on the relevant points and
                // calculate the area-average value on the face
                forAll(faces[facei], facePointi)
                {
                    const label pointi = faces[facei][facePointi];
                    pointFs[pointi] =
                        (pointXs[pointi] - cutXs[cutj - 1])
                       /(cutXs[cutj] - cutXs[cutj - 1]);
                }
                const scalar faceF =
                    cutPoly::faceAreaAverage
                    (
                        faces[facei],
                        points,
                        pointFs
                    ).second();

                // Add the whole face's contribution
                dynWeights.last().value += faceF*a;

                // Cut off anything before the left point
                if (faceMinXs[facei] < cutXs[cutj - 1])
                {
                    dynWeights.last().value -=
                        cutPoly::faceCutAreaIntegral
                        (
                            faces[facei],
                            faceAreas[facei],
                            faceF,
                            cutPoly::faceCuts
                            (
                                faces[facei],
                                pointXs,
                                cutXs[cutj - 1]
                            ),
                            points,
                            pointFs,
                            pointXs,
                            cutXs[cutj - 1],
                            true
                        ).second() & faceNormals[facei];
                }

                // Cut off anything after the middle point
                if (faceMaxXs[facei] > cutXs[cutj])
                {
                    dynWeights.last().value -=
                        cutPoly::faceCutAreaIntegral
                        (
                            faces[facei],
                            faceAreas[facei],
                            faceF,
                            cutPoly::faceCuts
                            (
                                faces[facei],
                                pointXs,
                                cutXs[cutj]
                            ),
                            points,
                            pointFs,
                            pointXs,
                            cutXs[cutj],
                            false
                        ).second() & faceNormals[facei];
                }
            }

            // Right interval
            if (cutj < cutXs.size() - 1 && faceMaxXs[facei] > cutXs[cutj])
            {
                // Update the basis function on the relevant points and
                // calculate the area-average value on the face
                forAll(faces[facei], facePointi)
                {
                    const label pointi = faces[facei][facePointi];
                    pointFs[pointi] =
                        (cutXs[cutj + 1] - pointXs[pointi])
                       /(cutXs[cutj + 1] - cutXs[cutj]);
                }
                const scalar faceF =
                    cutPoly::faceAreaAverage
                    (
                        faces[facei],
                        points,
                        pointFs
                    ).second();

                // Add the whole face's contribution
                dynWeights.last().value += faceF*a;

                // Cut off anything before the middle point
                if (faceMinXs[facei] < cutXs[cutj])
                {
                    dynWeights.last().value -=
                        cutPoly::faceCutAreaIntegral
                        (
                            faces[facei],
                            faceAreas[facei],
                            faceF,
                            cutPoly::faceCuts
                            (
                                faces[facei],
                                pointXs,
                                cutXs[cutj]
                            ),
                            points,
                            pointFs,
                            pointXs,
                            cutXs[cutj],
                            true
                        ).second() & faceNormals[facei];
                }

                // Cut off anything after the right point
                if (faceMaxXs[facei] > cutXs[cutj + 1])
                {
                    dynWeights.last().value -=
                        cutPoly::faceCutAreaIntegral
                        (
                            faces[facei],
                            faceAreas[facei],
                            faceF,
                            cutPoly::faceCuts
                            (
                                faces[facei],
                                pointXs,
                                cutXs[cutj + 1]
                            ),
                            points,
                            pointFs,
                            pointXs,
                            cutXs[cutj + 1],
                            false
                        ).second() & faceNormals[facei];
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
    const faceList& faces,
    const UList<vector>& faceAreas,
    const UList<vector>& faceNormals,
    const pointField& points,
    const scalarField& pointXs,
    const scalarField& faceMinXs,
    const scalarField& faceMaxXs,
    const labelList& faceMinOrder,
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
            faces,
            faceAreas,
            faceNormals,
            points,
            pointXs,
            faceMinXs,
            faceMaxXs,
            faceMinOrder,
            cutXs,
            normalise
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchCutPlot
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::List<Foam::patchCutPlot::weight> Foam::patchCutPlot::calcWeights
(
    const faceList& faces,
    const UList<vector>& faceAreas,
    const UList<vector>& faceNormals,
    const pointField& points,
    const scalarField& pointXs,
    const scalarField& cutXs,
    const bool interpolate,
    const bool normalise
)
{
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

    // Calculate and return the weights
    return
        calcWeights
        (
            faces,
            faceAreas,
            faceNormals,
            points,
            pointXs,
            faceMinXs,
            faceMaxXs,
            faceMinOrder,
            cutXs,
            interpolate,
            normalise
        );
}


Foam::tmp<Foam::scalarField> Foam::patchCutPlot::calcCutXs
(
    const faceList& faces,
    const UList<vector>& faceAreas,
    const UList<vector>& faceNormals,
    const pointField& points,
    const scalarField& pointXs,
    const bool interpolate,
    const label nCuts,
    const label nIter,
    const bool debug,
    const word& functionName,
    const polyMesh& functionMesh,
    const setWriter& functionFormatter
)
{
    const label nIntervals = interpolate ? nCuts : nCuts - 1;

    const scalarField faceMagAreas(mag(faceAreas));

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
        const List<patchCutPlot::weight> weights =
            patchCutPlot::calcWeights
            (
                faces,
                faceAreas,
                faceNormals,
                points,
                pointXs,
                faceMinXs,
                faceMaxXs,
                faceMinOrder,
                cutXs,
                interpolate,
                false
            );
        const scalarField intervalCounts
        (
            cutPlot::applyWeights(nIntervals, weights, faceMagAreas)
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
            functionMesh.time().globalPath()
           /functionObjects::writeFile::outputPrefix
           /(
               functionMesh.name() != polyMesh::defaultRegion
             ? functionMesh.name()
             : word()
            )
           /functionName
           /functionMesh.time().name();

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

    // Finally, calculate and return the actual normalised weights
    return tcutXs;
}


void Foam::patchCutPlot::writeLayers
(
    const SubList<face>& faces,
    const List<patchCutPlot::weight>& weights,
    const word& functionName,
    const fvMesh& functionMesh
)
{
    volTensorField layers
    (
        IOobject
        (
            functionName + ":layers",
            functionMesh.time().name(),
            functionMesh
        ),
        functionMesh,
        dimensionedTensor(dimless, tensor::uniform(-1))
    );

    label meshFacei0 = -1;
    forAll(functionMesh.faces(), meshFacei)
    {
        if (&functionMesh.faces()[meshFacei] == &faces[0])
        {
            meshFacei0 = meshFacei;
            break;
        }
    }

    forAll(weights, weighti)
    {
        const patchCutPlot::weight& w = weights[weighti];

        const label facei = meshFacei0 + w.elementi;
        const label patchi =
            functionMesh.boundaryMesh().patchIndices()
            [w.elementi + meshFacei0 - functionMesh.nInternalFaces()];
        const label patchFacei =
            facei - functionMesh.boundaryMesh()[patchi].start();

        layers.boundaryFieldRef()[patchi][patchFacei] = tensor::zero;
    }

    forAll(weights, weighti)
    {
        const patchCutPlot::weight& w = weights[weighti];

        const label facei = meshFacei0 + w.elementi;
        const label patchi =
            functionMesh.boundaryMesh().patchIndices()
            [w.elementi + meshFacei0 - functionMesh.nInternalFaces()];
        const label patchFacei =
            facei - functionMesh.boundaryMesh()[patchi].start();

        const direction i = w.cuti % tensor::nComponents;

        layers.boundaryFieldRef()[patchi][patchFacei][i] =
            w.value/functionMesh.magFaceAreas()[facei];
    }

    Info<< functionName << ": Writing " << layers.name() << endl;

    layers.write();
}


void Foam::patchCutPlot::writeLayers
(
    const UIndirectList<face>& faces,
    const List<patchCutPlot::weight>& weights,
    const word& functionName,
    const fvMesh& functionMesh
)
{
    volTensorField layers
    (
        IOobject
        (
            functionName + ":layers",
            functionMesh.time().name(),
            functionMesh
        ),
        functionMesh,
        dimensionedTensor(dimless, tensor::uniform(-1))
    );

    forAll(weights, weighti)
    {
        const patchCutPlot::weight& w = weights[weighti];

        const label facei = faces.addressing()[w.elementi];
        const label patchi =
            functionMesh.boundaryMesh().patchIndices()
            [w.elementi - functionMesh.nInternalFaces()];
        const label patchFacei =
            facei - functionMesh.boundaryMesh()[patchi].start();

        layers.boundaryFieldRef()[patchi][patchFacei] = tensor::zero;
    }

    forAll(weights, weighti)
    {
        const patchCutPlot::weight& w = weights[weighti];

        const label facei = faces.addressing()[w.elementi];
        const label patchi =
            functionMesh.boundaryMesh().patchIndices()
            [w.elementi - functionMesh.nInternalFaces()];
        const label patchFacei =
            facei - functionMesh.boundaryMesh()[patchi].start();

        const direction i = w.cuti % tensor::nComponents;

        layers.boundaryFieldRef()[patchi][patchFacei][i] =
            w.value/functionMesh.magFaceAreas()[facei];
    }

    Info<< functionName << ": Writing " << layers.name() << endl;

    layers.write();
}


// ************************************************************************* //
