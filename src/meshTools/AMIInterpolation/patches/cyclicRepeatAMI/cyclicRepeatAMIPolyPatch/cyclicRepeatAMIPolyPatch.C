/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "cyclicRepeatAMIPolyPatch.H"
#include "SubField.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicRepeatAMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicRepeatAMIPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, cyclicRepeatAMIPolyPatch, dictionary);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template <class Type>
class keepIfTrueOp
{
    public:

        Tuple2<bool, Type> operator()
        (
            const Tuple2<bool, Type>& x,
            const Tuple2<bool, Type>& y
        ) const
        {
            if (x.first())
            {
                return x;
            }
            else if (y.first())
            {
                return y;
            }
            else
            {
                return Tuple2<bool, Type>(false, Type());
            }
        }
};

}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::cyclicRepeatAMIPolyPatch::resetAMI() const
{
    if (!owner())
    {
        return;
    }

    Info<< indent << typeName <<" : Creating addressing and weights between "
        << this->size() << " source faces and " << neighbPatch().size()
        << " target faces" << endl;

    // Get the transform associated with the transform patch
    vectorTensorTransform t;
    {
        const coupledPolyPatch& transformPatch = this->transformPatch();

        if (transformPatch.name() != neighbPatch().transformPatch().name())
        {
            FatalErrorInFunction
                << "Transform patch " << transformPatch.name() << " for "
                << typeName << " patch " << name() << " is not the same as for "
                << "the neighbour patch " << neighbPatch().name() << ". "
                << "This is not allowed." << exit(FatalError);
        }

        if
        (
            transformPatch.separation().size() > 1
         || transformPatch.forwardT().size() > 1
        )
        {
            FatalErrorInFunction
                << "Transform patch " << transformPatch.name() << " for "
                << typeName << " patch " << name() << " has a non-uniform "
                << "transformation. This is not allowed."
                << exit(FatalError);
        }

        Tuple2<bool, vectorTensorTransform> bt
        (
            transformPatch.size(),
            vectorTensorTransform
            (
                transformPatch.separation().size() > 0
              ? transformPatch.separation()[0]
              : vector::zero,
                transformPatch.forwardT().size() > 0
              ? transformPatch.forwardT()[0]
              : tensor::zero,
                transformPatch.forwardT().size() > 0
            )
        );

        reduce(bt, keepIfTrueOp<vectorTensorTransform>());

        t = bt.second();
    }
    const vectorTensorTransform tInv(inv(t));

    // Work out the number of repetitions of the transform that separate this
    // patch from its neighbour
    label n = 0;
    {
        const scalarField thisMagAreas(mag(this->faceAreas()));
        const scalarField nbrMagAreas(mag(neighbPatch().faceAreas()));

        vector thisCentre =
            gSum(this->faceCentres()*thisMagAreas)/gSum(thisMagAreas);
        vector nbrCentre =
            gSum(neighbPatch().faceCentres()*nbrMagAreas)/gSum(nbrMagAreas);

        scalar dLeft = mag(t.transformPosition(thisCentre) - nbrCentre);
        scalar d = mag(thisCentre - nbrCentre);
        scalar dRight = mag(tInv.transformPosition(thisCentre) - nbrCentre);

        while (dLeft < d)
        {
            thisCentre = t.transformPosition(thisCentre);

            dRight = d;
            d = dLeft;
            dLeft = mag(t.transformPosition(thisCentre) - nbrCentre);

            ++ n;
        }

        while (dRight < d)
        {
            thisCentre = tInv.transformPosition(thisCentre);

            dLeft = d;
            d = dRight;
            dRight = mag(tInv.transformPosition(thisCentre) - nbrCentre);

            -- n;
        }
    }

    // Generate the full transformations
    vectorTensorTransform TLeft(t), T(vectorTensorTransform::I), TRight(tInv);
    if (n > 0)
    {
        for (label i = 0; i < n - 1; ++ i)
        {
            T = t & T;
        }

        TLeft = T;
        T = t & T;
        TRight = t & T;
    }
    if (n < 0)
    {
        for (label i = 0; i > n + 1; -- i)
        {
            T = tInv & T;
        }

        TRight = T;
        T = tInv & T;
        TLeft = tInv & T;
    }

    // Create copies of this patch and the neighbour patch's points
    pointField thisPoints(localPoints());
    const pointField nbrPoints(neighbPatch().localPoints());

    // Create primitive patches
    primitivePatch thisPatch
    (
        SubList<face>(localFaces(), size()),
        thisPoints
    );
    primitivePatch nbrPatch
    (
        SubList<face>(neighbPatch().localFaces(), neighbPatch().size()),
        nbrPoints
    );

    // Do the three bounding AMI interpolations
    thisPoints = TLeft.transformPosition(localPoints());
    thisPatch.movePoints(thisPoints);
    autoPtr<AMIPatchToPatchInterpolation> AMILeftPtr
    (
        new AMIPatchToPatchInterpolation
        (
            thisPatch,
            nbrPatch,
            faceAreaIntersect::tmMesh,
            false,
            AMIPatchToPatchInterpolation::imPartialFaceAreaWeight,
            AMILowWeightCorrection_,
            AMIReverse_,
            false
        )
    );
    const scalar sLeft =
        gSum(AMILeftPtr->srcWeightsSum()*AMILeftPtr->srcMagSf())
       /gSum(AMILeftPtr->srcMagSf());

    thisPoints = T.transformPosition(localPoints());
    thisPatch.movePoints(thisPoints);
    autoPtr<AMIPatchToPatchInterpolation> AMIPtr
    (
        new AMIPatchToPatchInterpolation
        (
            thisPatch,
            nbrPatch,
            faceAreaIntersect::tmMesh,
            false,
            AMIPatchToPatchInterpolation::imPartialFaceAreaWeight,
            AMILowWeightCorrection_,
            AMIReverse_,
            false
        )
    );
    const scalar s =
        gSum(AMIPtr->srcWeightsSum()*AMIPtr->srcMagSf())
       /gSum(AMIPtr->srcMagSf());

    thisPoints = TRight.transformPosition(localPoints());
    thisPatch.movePoints(thisPoints);
    autoPtr<AMIPatchToPatchInterpolation> AMIRightPtr
    (
        new AMIPatchToPatchInterpolation
        (
            thisPatch,
            nbrPatch,
            faceAreaIntersect::tmMesh,
            false,
            AMIPatchToPatchInterpolation::imPartialFaceAreaWeight,
            AMILowWeightCorrection_,
            AMIReverse_,
            false
        )
    );
    const scalar sRight =
        gSum(AMIRightPtr->srcWeightsSum()*AMIRightPtr->srcMagSf())
       /gSum(AMIRightPtr->srcMagSf());

    Info<< typeName << ": number of transforms = " << n << endl
        << typeName << ": left/centre/right sum(weights) = "
        << sLeft << ", " << s << ", " << sRight << endl;

    // Set the AMI interpolators and transforms using the centre and the most
    // overlapping of the left and right sides
    AMIs_.resize(2);
    AMIs_.set(0, sLeft > sRight ? AMILeftPtr.ptr() : AMIPtr.ptr());
    AMIs_.set(1, sLeft > sRight ? AMIPtr.ptr() : AMIRightPtr.ptr());

    AMITransforms_.resize(2);
    AMITransforms_[0] = sLeft > sRight ? TLeft : T;
    AMITransforms_[1] = sLeft > sRight ? T : TRight;

    // Sum and normalise the two AMI interpolators
    AMIPatchToPatchInterpolation::sumWeights(AMIs_);
    AMIPatchToPatchInterpolation::reportSumWeights(AMIs_[0]);
    AMIPatchToPatchInterpolation::normaliseWeights(AMIs_);
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::cyclicRepeatAMIPolyPatch::cyclicRepeatAMIPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const transformType transform
)
:
    cyclicAMIPolyPatch
    (
        name,
        size,
        start,
        index,
        bm,
        patchType,
        transform,
        false,
        AMIPatchToPatchInterpolation::imFaceAreaWeight
    ),
    transformPatchName_(word::null),
    transformPatchID_(-1)
{
    // Transform patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicRepeatAMIPolyPatch::cyclicRepeatAMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicAMIPolyPatch
    (
        name,
        dict,
        index,
        bm,
        patchType,
        false,
        AMIPatchToPatchInterpolation::imFaceAreaWeight
    ),
    transformPatchName_(dict.lookup("transformPatch")),
    transformPatchID_(-1)
{
    // Transform patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicRepeatAMIPolyPatch::cyclicRepeatAMIPolyPatch
(
    const cyclicRepeatAMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    cyclicAMIPolyPatch(pp, bm),
    transformPatchName_(pp.transformPatchName_),
    transformPatchID_(-1)
{
    // Transform patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicRepeatAMIPolyPatch::cyclicRepeatAMIPolyPatch
(
    const cyclicRepeatAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName
)
:
    cyclicAMIPolyPatch(pp, bm, index, newSize, newStart, nbrPatchName),
    transformPatchName_(pp.transformPatchName_),
    transformPatchID_(-1)
{
    // Transform patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicRepeatAMIPolyPatch::cyclicRepeatAMIPolyPatch
(
    const cyclicRepeatAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    cyclicAMIPolyPatch(pp, bm, index, mapAddressing, newStart),
    transformPatchName_(pp.transformPatchName_),
    transformPatchID_(-1)
{
    // Transform patch might not be valid yet so cannot determine
    // associated patchID
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicRepeatAMIPolyPatch::~cyclicRepeatAMIPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::cyclicRepeatAMIPolyPatch&
Foam::cyclicRepeatAMIPolyPatch::neighbPatch() const
{
    const polyPatch& pp = this->boundaryMesh()[neighbPatchID()];

    return refCast<const cyclicRepeatAMIPolyPatch>(pp);
}


Foam::label Foam::cyclicRepeatAMIPolyPatch::transformPatchID() const
{
    if (transformPatchID_ == -1)
    {
        transformPatchID_ =
            this->boundaryMesh().findPatchID(transformPatchName());

        if (transformPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Illegal transformPatch name " << transformPatchName()
                << nl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }
    }

    return transformPatchID_;
}


const Foam::coupledPolyPatch&
Foam::cyclicRepeatAMIPolyPatch::transformPatch() const
{
    const polyPatch& pp = this->boundaryMesh()[transformPatchID()];

    return refCast<const coupledPolyPatch>(pp);
}


const Foam::scalarField& Foam::cyclicRepeatAMIPolyPatch::weightsSum() const
{
    // The two AMI-interpolation classes have their weights summed together, so
    // both should contain the same weights sum field. We can, therefore
    // delegate to the base class and just return the weights sum of the first.

    return cyclicAMIPolyPatch::weightsSum();
}


const Foam::scalarField&
Foam::cyclicRepeatAMIPolyPatch::neighbWeightsSum() const
{
    // See above.

    return cyclicAMIPolyPatch::neighbWeightsSum();
}


void Foam::cyclicRepeatAMIPolyPatch::write(Ostream& os) const
{
    cyclicAMIPolyPatch::write(os);

    os.writeKeyword("transformPatch") << transformPatchName_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //
