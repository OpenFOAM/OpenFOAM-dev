/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "globalIndexAndTransform.H"
#include "cyclicPolyPatch.H"
#include "DynamicField.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * Private Static Data Members * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(globalIndexAndTransform, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::globalIndexAndTransform::matchTransform
(
    const List<vectorTensorTransform>& refTransforms,
    label& matchedRefTransformI,
    const vectorTensorTransform& testTransform,
    scalar tolerance,
    bool checkBothSigns
) const
{
    matchedRefTransformI = -1;

    forAll(refTransforms, i)
    {
        const vectorTensorTransform& refTransform = refTransforms[i];

        scalar maxVectorMag = sqrt
        (
            max(magSqr(testTransform.t()), magSqr(refTransform.t()))
        );

        // Test the difference between vector parts to see if it is
        // less than tolerance times the larger vector part magnitude.

        scalar vectorDiff =
            mag(refTransform.t() - testTransform.t())
           /(maxVectorMag + vSmall)
           /tolerance;

        // Test the difference between tensor parts to see if it is
        // less than the tolerance.  sqrt(3.0) factor used to scale
        // differences as this is magnitude of a rotation tensor.  If
        // neither transform has a rotation, then the test is not
        // necessary.

        scalar tensorDiff = 0;

        if (refTransform.hasR() || testTransform.hasR())
        {
            tensorDiff =
                mag(refTransform.R() - testTransform.R())
               /sqrt(3.0)
               /tolerance;
        }

        // ...Diff result is < 1 if the test part matches the ref part
        // within tolerance

        if (vectorDiff < 1 && tensorDiff < 1)
        {
            matchedRefTransformI = i;

            return +1;
        }

        if (checkBothSigns)
        {
            // Test the inverse transform differences too

            vectorDiff =
                mag(refTransform.t() + testTransform.t())
               /(maxVectorMag + vSmall)
               /tolerance;

            tensorDiff = 0;

            if (refTransform.hasR() || testTransform.hasR())
            {
                tensorDiff =
                    mag(refTransform.R() - testTransform.R().T())
                   /sqrt(3.0)
                   /tolerance;
            }

            if (vectorDiff < 1 && tensorDiff < 1)
            {
                matchedRefTransformI = i;

                return -1;
            }
        }
    }

    return 0;
}


void Foam::globalIndexAndTransform::determineTransforms()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    DynamicList<vectorTensorTransform> localTransforms;
    DynamicField<scalar> localTols;

    label dummyMatch = -1;

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        // Note: special check for unordered cyclics. These are in fact
        // transform bcs and should probably be split off.
        // Note: We don't want to be finding transforms for patches marked as
        // coincident full match. These should have no transform by definition.
        if
        (
            isA<coupledPolyPatch>(pp)
        && !(
                isA<cyclicPolyPatch>(pp)
             && refCast<const cyclicPolyPatch>(pp).transform()
             == cyclicPolyPatch::NOORDERING
            )
        && !(
                refCast<const coupledPolyPatch>(pp).transform()
             == coupledPolyPatch::COINCIDENTFULLMATCH
            )
        )
        {
            const coupledPolyPatch& cpp = refCast<const coupledPolyPatch>(pp);

            if (cpp.separated())
            {
                const vectorField& sepVecs = cpp.separation();

                forAll(sepVecs, sVI)
                {
                    const vector& sepVec = sepVecs[sVI];

                    if (mag(sepVec) > small)
                    {
                        vectorTensorTransform transform(sepVec);

                        if
                        (
                            matchTransform
                            (
                                localTransforms,
                                dummyMatch,
                                transform,
                                cpp.matchTolerance(),
                                false
                            ) == 0
                        )
                        {
                            localTransforms.append(transform);
                            localTols.append(cpp.matchTolerance());
                        }
                    }
                }
            }
            else if (!cpp.parallel())
            {
                const tensorField& transTensors = cpp.reverseT();

                forAll(transTensors, tTI)
                {
                    const tensor& transT = transTensors[tTI];

                    if (mag(transT - I) > small)
                    {
                        vectorTensorTransform transform(transT);

                        if
                        (
                            matchTransform
                            (
                                localTransforms,
                                dummyMatch,
                                transform,
                                cpp.matchTolerance(),
                                false
                            ) == 0
                        )
                        {
                            localTransforms.append(transform);
                            localTols.append(cpp.matchTolerance());
                        }
                    }
                }
            }
        }
    }


    // Collect transforms on master
    List<List<vectorTensorTransform>> allTransforms(Pstream::nProcs());
    allTransforms[Pstream::myProcNo()] = localTransforms;
    Pstream::gatherList(allTransforms);

    // Collect matching tolerance on master
    List<scalarField> allTols(Pstream::nProcs());
    allTols[Pstream::myProcNo()] = localTols;
    Pstream::gatherList(allTols);

    if (Pstream::master())
    {
        localTransforms.clear();

        forAll(allTransforms, proci)
        {
            const List<vectorTensorTransform>& procTransVecs =
                allTransforms[proci];

            forAll(procTransVecs, pSVI)
            {
                const vectorTensorTransform& transform = procTransVecs[pSVI];

                if (mag(transform.t()) > small || transform.hasR())
                {
                    if
                    (
                        matchTransform
                        (
                            localTransforms,
                            dummyMatch,
                            transform,
                            allTols[proci][pSVI],
                            true
                        ) == 0
                    )
                    {
                        localTransforms.append(transform);
                    }
                }
            }
        }
    }

    transforms_.transfer(localTransforms);
    Pstream::scatter(transforms_);
}


void Foam::globalIndexAndTransform::determineTransformPermutations()
{
    label nTransformPermutations = pow(label(3), transforms_.size());

    transformPermutations_.setSize(nTransformPermutations);

    forAll(transformPermutations_, tPI)
    {
        vectorTensorTransform transform;

        label transformIndex = tPI;

        // Invert the ternary index encoding using repeated division by
        // three

        forAll(transforms_, b)
        {
            const label w = (transformIndex % 3) - 1;

            transformIndex /= 3;

            if (w > 0)
            {
                transform &= transforms_[b];
            }
            else if (w < 0)
            {
                transform &= inv(transforms_[b]);
            }
        }

        transformPermutations_[tPI] = transform;
    }


    // Encode index for 0 sign
    labelList permutationIndices(nIndependentTransforms(), 0);
    nullTransformIndex_ = encodeTransformIndex(permutationIndices);
}


void Foam::globalIndexAndTransform::determinePatchTransformSign()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    patchTransformSign_.setSize(patches.size(), labelPair(-1, 0));

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        // Note: special check for unordered cyclics. These are in fact
        // transform bcs and should probably be split off.
        // Note: We don't want to be finding transforms for patches marked as
        // coincident full match. These should have no transform by definition.
        if
        (
            isA<coupledPolyPatch>(pp)
        && !(
                isA<cyclicPolyPatch>(pp)
             && refCast<const cyclicPolyPatch>(pp).transform()
             == cyclicPolyPatch::NOORDERING
            )
        && !(
                refCast<const coupledPolyPatch>(pp).transform()
             == coupledPolyPatch::COINCIDENTFULLMATCH
            )
        )
        {
            const coupledPolyPatch& cpp = refCast<const coupledPolyPatch>(pp);

            if (cpp.separated())
            {
                const vectorField& sepVecs = cpp.separation();

                // This loop is implicitly expecting only a single
                // value for separation()
                forAll(sepVecs, sVI)
                {
                    const vector& sepVec = sepVecs[sVI];

                    if (mag(sepVec) > small)
                    {
                        vectorTensorTransform t(sepVec);

                        label matchTransI;
                        label sign = matchTransform
                        (
                            transforms_,
                            matchTransI,
                            t,
                            cpp.matchTolerance(),
                            true
                        );
                        patchTransformSign_[patchi] =
                            labelPair(matchTransI, sign);
                    }
                }

            }
            else if (!cpp.parallel())
            {
                const tensorField& transTensors = cpp.reverseT();

                // This loop is implicitly expecting only a single
                // value for reverseT()
                forAll(transTensors, tTI)
                {
                    const tensor& transT = transTensors[tTI];

                    if (mag(transT - I) > small)
                    {
                        vectorTensorTransform t(transT);

                        label matchTransI;
                        label sign = matchTransform
                        (
                            transforms_,
                            matchTransI,
                            t,
                            cpp.matchTolerance(),
                            true
                        );

                        patchTransformSign_[patchi] =
                            labelPair(matchTransI, sign);
                    }
                }
            }
        }
    }
}


bool Foam::globalIndexAndTransform::uniqueTransform
(
    const point& pt,
    labelPairList& trafos,
    const label patchi,
    const labelPair& patchTrafo
) const
{
    if (findIndex(trafos, patchTrafo) == -1)
    {
        // New transform. Check if already have 3
        if (trafos.size() == 3)
        {
            if (patchi > -1)
            {
                WarningInFunction
                    << "Point " << pt
                    << " is on patch " << mesh_.boundaryMesh()[patchi].name();
            }
            else
            {
                WarningInFunction
                    << "Point " << pt << " is on a coupled patch";
            }
            Warning
                << " with transformation " << patchTrafo
                << " but also on 3 other patches with transforms "
                << trafos << nl
                << "This is not a space filling tiling and might"
                << " indicate a setup problem and give problems"
                << " for e.g. lagrangian tracking or interpolation" << endl;

            // Already warned so no need to extend more
            trafos.clear();
            return false;
        }

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalIndexAndTransform::globalIndexAndTransform(const polyMesh& mesh)
:
    mesh_(mesh),
    transforms_(),
    transformPermutations_(),
    patchTransformSign_()
{
    determineTransforms();

    determineTransformPermutations();

    determinePatchTransformSign();

    if (debug && transforms_.size() > 0)
    {
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        Info<< "Determined global transforms :" << endl;
        Info<< "\t\ttranslation\trotation" << endl;
        forAll(transforms_, i)
        {
            Info<< '\t' << i << '\t';
            const vectorTensorTransform& trafo = transforms_[i];
            if (trafo.hasR())
            {
                 Info<< trafo.t() << '\t' << trafo.R();
            }
            else
            {
                 Info<< trafo.t() << '\t' << "---";
            }
            Info<< endl;
        }
        Info<< endl;


        Info<< "\tpatch\ttransform\tsign" << endl;
        forAll(patchTransformSign_, patchi)
        {
            if (patchTransformSign_[patchi].first() != -1)
            {
                Info<< '\t' << patches[patchi].name()
                    << '\t' << patchTransformSign_[patchi].first()
                    << '\t' << patchTransformSign_[patchi].second()
                    << endl;
            }
        }
        Info<< endl;


        Info<< "Permutations of transformations:" << endl
            << "\t\ttranslation\trotation" << endl;
        forAll(transformPermutations_, i)
        {
            Info<< '\t' << i << '\t';
            const vectorTensorTransform& trafo = transformPermutations_[i];
            if (trafo.hasR())
            {
                 Info<< trafo.t() << '\t' << trafo.R();
            }
            else
            {
                 Info<< trafo.t() << '\t' << "---";
            }
            Info<< endl;
        }
        Info<< "nullTransformIndex:" << nullTransformIndex() << endl
            << endl;
    }


    if (transforms_.size() > 0)
    {
        // Check that the transforms are space filling : any point
        // can only use up to three transforms

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();


        // 1. Collect transform&sign per point and do local check

        List<labelPairList> pointToTrafos(mesh_.nPoints());

        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            const labelPair& transSign = patchTransformSign_[patchi];

            if (transSign.first() > -1)
            {
                const labelList& mp = pp.meshPoints();
                forAll(mp, i)
                {
                    labelPairList& trafos = pointToTrafos[mp[i]];

                    bool newTransform = uniqueTransform
                    (
                        mesh_.points()[mp[i]],
                        trafos,
                        patchi,
                        transSign
                    );

                    if (newTransform)
                    {
                        trafos.append(transSign);
                    }
                }
            }
        }


        // Synchronise across collocated (= untransformed) points
        // TBD: there is a big problem in that globalMeshData uses
        //      globalIndexAndTransform. Triggers recursion.
        if (false)
        {
            const globalMeshData& gmd = mesh_.globalData();
            const indirectPrimitivePatch& cpp = gmd.coupledPatch();
            const labelList& meshPoints = cpp.meshPoints();
            const mapDistribute& slavesMap = gmd.globalCoPointSlavesMap();
            const labelListList& slaves = gmd.globalCoPointSlaves();

            List<labelPairList> elems(slavesMap.constructSize());
            forAll(meshPoints, i)
            {
                elems[i] = pointToTrafos[meshPoints[i]];
            }

            // Pull slave data onto master. No need to update transformed slots.
            slavesMap.distribute(elems, false);

            // Combine master data with slave data
            forAll(slaves, i)
            {
                labelPairList& trafos = elems[i];

                const labelList& slavePoints = slaves[i];

                // Combine master with untransformed slave data
                forAll(slavePoints, j)
                {
                    const labelPairList& slaveTrafos = elems[slavePoints[j]];

                    forAll(slaveTrafos, slaveI)
                    {
                        bool newTransform = uniqueTransform
                        (
                            mesh_.points()[meshPoints[i]],
                            trafos,
                            -1,
                            slaveTrafos[slaveI]
                        );

                        if (newTransform)
                        {
                            trafos.append(slaveTrafos[slaveI]);
                        }
                    }
                }
            }
        }
    }
}


// ************************************************************************* //
