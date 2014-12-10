/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

// * * * * * * * * * * * * Private Static Data Members * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(globalIndexAndTransform, 0);
const label globalIndexAndTransform::base_ = 32;
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
           /(maxVectorMag + VSMALL)
           /tolerance;

        // Test the difference between tensor parts to see if it is
        // less than the tolerance.  sqrt(3.0) factor used to scale
        // differnces as this is magnitude of a rotation tensor.  If
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
               /(maxVectorMag + VSMALL)
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

    transforms_ = List<vectorTensorTransform>(6);
    scalarField maxTol(6);

    label nextTrans = 0;

    label dummyMatch = -1;

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        // Note: special check for unordered cyclics. These are in fact
        // transform bcs and should probably be split off.
        if
        (
            isA<coupledPolyPatch>(pp)
        && !(
                isA<cyclicPolyPatch>(pp)
             && (
                    refCast<const cyclicPolyPatch>(pp).transform()
                 == cyclicPolyPatch::NOORDERING
                )
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

                    if (mag(sepVec) > SMALL)
                    {
                        vectorTensorTransform transform(sepVec);

                        if
                        (
                            matchTransform
                            (
                                transforms_,
                                dummyMatch,
                                transform,
                                cpp.matchTolerance(),
                                false
                            ) == 0
                        )
                        {
                            if (nextTrans == 6)
                            {
                                FatalErrorIn
                                (
                                     "void Foam::globalIndexAndTransform::"
                                     "determineTransforms()"
                                )   << "More than six unsigned transforms"
                                    << " detected:" << nl << transforms_
                                    << exit(FatalError);
                            }
                            transforms_[nextTrans] = transform;
                            maxTol[nextTrans++] = cpp.matchTolerance();
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

                    if (mag(transT - I) > SMALL)
                    {
                        vectorTensorTransform transform(transT);

                        if
                        (
                            matchTransform
                            (
                                transforms_,
                                dummyMatch,
                                transform,
                                cpp.matchTolerance(),
                                false
                            ) == 0
                        )
                        {
                            if (nextTrans == 6)
                            {
                                FatalErrorIn
                                (
                                    "void Foam::globalIndexAndTransform::"
                                    "determineTransforms()"
                                )   << "More than six unsigned transforms"
                                    << " detected:" << nl << transforms_
                                    << exit(FatalError);
                            }
                            transforms_[nextTrans] = transform;
                            maxTol[nextTrans++] = cpp.matchTolerance();
                        }
                    }
                }
            }
        }
    }


    // Collect transforms on master

    List<List<vectorTensorTransform> > allTransforms(Pstream::nProcs());
    allTransforms[Pstream::myProcNo()] = transforms_;
    Pstream::gatherList(allTransforms);

    // Collect matching tolerance on master
    List<scalarField> allTols(Pstream::nProcs());
    allTols[Pstream::myProcNo()] = maxTol;
    Pstream::gatherList(allTols);

    if (Pstream::master())
    {
        transforms_ = List<vectorTensorTransform>(3);

        label nextTrans = 0;

        forAll(allTransforms, procI)
        {
            const List<vectorTensorTransform>& procTransVecs =
                allTransforms[procI];

            forAll(procTransVecs, pSVI)
            {
                const vectorTensorTransform& transform = procTransVecs[pSVI];

                if (mag(transform.t()) > SMALL || transform.hasR())
                {
                    if
                    (
                        matchTransform
                        (
                            transforms_,
                            dummyMatch,
                            transform,
                            allTols[procI][pSVI],
                            true
                        ) ==  0
                    )
                    {
                        transforms_[nextTrans++] = transform;
                    }

                    if (nextTrans > 3)
                    {
                        FatalErrorIn
                        (
                            "void Foam::globalIndexAndTransform::"
                            "determineTransforms()"
                        )
                            << "More than three independent basic "
                            << "transforms detected:" << nl
                            << allTransforms
                            << transforms_
                            << exit(FatalError);
                    }
                }
            }
        }

        transforms_.setSize(nextTrans);
    }

    Pstream::scatter(transforms_);

    if (transforms_.size() > 3)
    {
        WarningIn
        (
            "void globalIndexAndTransform::determineTransforms()"
        )   << "More than three independent basic "
            << "transforms detected:" << nl
            << transforms_ << nl
            << "This is not a space filling tiling and will probably"
            << " give problems for e.g. lagrangian tracking or interpolation"
            << endl;
    }
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

    patchTransformSign_.setSize(patches.size(), Pair<label>(-1, 0));

    label matchTransI = -1;

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        // Pout<< nl << patchI << " " << pp.name() << endl;

        // Note: special check for unordered cyclics. These are in fact
        // transform bcs and should probably be split off.
        if
        (
            isA<coupledPolyPatch>(pp)
        && !(
                isA<cyclicPolyPatch>(pp)
             && (
                    refCast<const cyclicPolyPatch>(pp).transform()
                 == cyclicPolyPatch::NOORDERING
                )
            )
        )
        {
            const coupledPolyPatch& cpp =
            refCast<const coupledPolyPatch>(pp);

            if (cpp.separated())
            {
                const vectorField& sepVecs = cpp.separation();

                // Pout<< "sepVecs " << sepVecs << endl;

                // This loop is implicitly expecting only a single
                // value for separation()
                forAll(sepVecs, sVI)
                {
                    const vector& sepVec = sepVecs[sVI];

                    if (mag(sepVec) > SMALL)
                    {
                        vectorTensorTransform t(sepVec);

                        label sign = matchTransform
                        (
                            transforms_,
                            matchTransI,
                            t,
                            cpp.matchTolerance(),
                            true
                        );

                        // Pout<< sign << " " << matchTransI << endl;

                        // List<label> permutation(transforms_.size(), 0);

                        // permutation[matchTransI] = sign;

                        // Pout<< encodeTransformIndex(permutation) << nl
                        //     << transformPermutations_
                        //        [
                        //            encodeTransformIndex(permutation)
                        //        ]
                        //     << endl;

                        patchTransformSign_[patchI] =
                            Pair<label>(matchTransI, sign);
                    }
                }

            }
            else if (!cpp.parallel())
            {
                const tensorField& transTensors = cpp.reverseT();

                // Pout<< "transTensors " << transTensors << endl;

                // This loop is implicitly expecting only a single
                // value for reverseT()
                forAll(transTensors, tTI)
                {
                    const tensor& transT = transTensors[tTI];

                    if (mag(transT - I) > SMALL)
                    {
                        vectorTensorTransform t(transT);

                        label sign = matchTransform
                        (
                            transforms_,
                            matchTransI,
                            t,
                            cpp.matchTolerance(),
                            true
                        );

                        // Pout<< sign << " " << matchTransI << endl;

                        // List<label> permutation(transforms_.size(), 0);

                        // permutation[matchTransI] = sign;

                        // Pout<< encodeTransformIndex(permutation) << nl
                        //     << transformPermutations_
                        //        [
                        //            encodeTransformIndex(permutation)
                        //        ]
                        //     << endl;

                        patchTransformSign_[patchI] =
                            Pair<label>(matchTransI, sign);
                    }
                }
            }
        }
    }

    // Pout<< patchTransformSign_ << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalIndexAndTransform::globalIndexAndTransform
(
    const polyMesh& mesh
)
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
        forAll(patchTransformSign_, patchI)
        {
            if (patchTransformSign_[patchI].first() != -1)
            {
                Info<< '\t' << patches[patchI].name()
                    << '\t' << patchTransformSign_[patchI].first()
                    << '\t' << patchTransformSign_[patchI].second()
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::globalIndexAndTransform::~globalIndexAndTransform()
{}


// ************************************************************************* //
