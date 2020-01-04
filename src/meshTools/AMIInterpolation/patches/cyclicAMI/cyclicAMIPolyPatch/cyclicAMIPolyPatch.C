/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "cyclicAMIPolyPatch.H"
#include "SubField.H"
#include "Time.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicAMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicAMIPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, cyclicAMIPolyPatch, dictionary);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::vector Foam::cyclicAMIPolyPatch::findFaceNormalMaxRadius
(
    const pointField& faceCentres
) const
{
    // Determine a face furthest away from the axis

    const vectorField n((faceCentres - rotationCentre_) ^ rotationAxis_);

    const scalarField magRadSqr(magSqr(n));

    label facei = findMax(magRadSqr);

    if (debug)
    {
        Info<< "findFaceMaxRadius(const pointField&) : patch: " << name() << nl
            << "    rotFace  = " << facei << nl
            << "    point    = " << faceCentres[facei] << nl
            << "    distance = " << Foam::sqrt(magRadSqr[facei])
            << endl;
    }

    return n[facei];
}


void Foam::cyclicAMIPolyPatch::calcTransforms
(
    const primitivePatch& thisPatch,
    const pointField& thisPatchCtrs,
    const vectorField& thisPatchAreas,
    const pointField& nbrPatchCtrs,
    const vectorField& nbrPatchAreas
)
{
    if (transformType() != nbrPatch().transformType())
    {
        FatalErrorInFunction
            << "Patch " << name()
            << " has transform type " << transformTypeNames[transformType()]
            << ", neighbour patch " << nbrPatchName()
            << " has transform type "
            << nbrPatch().transformTypeNames[nbrPatch().transformType()]
            << exit(FatalError);
    }


    // Calculate transformation tensors

    switch (transformType())
    {
        case ROTATIONAL:
        {
            tensor revT = Zero;

            if (rotationAngleDefined_)
            {
                const tensor T(rotationAxis_*rotationAxis_);

                const tensor S
                (
                    0, -rotationAxis_.z(), rotationAxis_.y(),
                    rotationAxis_.z(), 0, -rotationAxis_.x(),
                    -rotationAxis_.y(), rotationAxis_.x(), 0
                );

                const tensor revTPos
                (
                    T
                  + cos(rotationAngle_)*(tensor::I - T)
                  + sin(rotationAngle_)*S
                );

                const tensor revTNeg
                (
                    T
                  + cos(-rotationAngle_)*(tensor::I - T)
                  + sin(-rotationAngle_)*S
                );

                // Check - assume correct angle when difference in face areas
                // is the smallest
                const vector transformedAreaPos = gSum(nbrPatchAreas & revTPos);
                const vector transformedAreaNeg = gSum(nbrPatchAreas & revTNeg);
                const vector area0 = gSum(thisPatchAreas);
                const scalar magArea0 = mag(area0) + rootVSmall;

                // Areas have opposite sign, so sum should be zero when correct
                // rotation applied
                const scalar errorPos = mag(transformedAreaPos + area0);
                const scalar errorNeg = mag(transformedAreaNeg + area0);

                const scalar normErrorPos = errorPos/magArea0;
                const scalar normErrorNeg = errorNeg/magArea0;

                if (errorPos > errorNeg && normErrorNeg < matchTolerance())
                {
                    revT = revTNeg;
                    rotationAngle_ *= -1;
                }
                else
                {
                    revT = revTPos;
                }

                const scalar areaError = min(normErrorPos, normErrorNeg);

                if (areaError > matchTolerance())
                {
                    WarningInFunction
                        << "Patch areas are not consistent within "
                        << 100*matchTolerance()
                        << " % indicating a possible error in the specified "
                        << "angle of rotation" << nl
                        << "    owner patch     : " << name() << nl
                        << "    neighbour patch : " << nbrPatch().name()
                        << nl
                        << "    angle           : "
                        << radToDeg(rotationAngle_) << " deg" << nl
                        << "    area error      : " << 100*areaError << " %"
                        << "    match tolerance : " <<  matchTolerance()
                        << endl;
                }

                if (debug)
                {
                    scalar theta = radToDeg(rotationAngle_);

                    Pout<< "cyclicAMIPolyPatch::calcTransforms: patch:"
                        << name()
                        << " Specified rotation:"
                        << " swept angle: " << theta << " [deg]"
                        << " reverse transform: " << revT
                        << endl;
                }
            }
            else
            {
                point n0 = Zero;
                point n1 = Zero;
                if (thisPatchCtrs.size())
                {
                    n0 = findFaceNormalMaxRadius(thisPatchCtrs);
                }
                if (nbrPatchCtrs.size())
                {
                    n1 = -findFaceNormalMaxRadius(nbrPatchCtrs);
                }

                reduce(n0, maxMagSqrOp<point>());
                reduce(n1, maxMagSqrOp<point>());

                n0 /= mag(n0) + vSmall;
                n1 /= mag(n1) + vSmall;

                // Extended tensor from two local coordinate systems calculated
                // using normal and rotation axis
                const tensor E0
                (
                    rotationAxis_,
                    (n0 ^ rotationAxis_),
                    n0
                );
                const tensor E1
                (
                    rotationAxis_,
                    (-n1 ^ rotationAxis_),
                    -n1
                );
                revT = E1.T() & E0;

                if (debug)
                {
                    scalar theta = radToDeg(acos(-(n0 & n1)));

                    Pout<< "cyclicAMIPolyPatch::calcTransforms: patch:"
                        << name()
                        << " Specified rotation:"
                        << " n0:" << n0 << " n1:" << n1
                        << " swept angle: " << theta << " [deg]"
                        << " reverse transform: " << revT
                        << endl;
                }
            }

            separation_ = Zero;

            transform_ = transformer(revT.T());

            break;
        }
        case TRANSLATIONAL:
        {
            if (debug)
            {
                Pout<< "cyclicAMIPolyPatch::calcTransforms : patch:" << name()
                    << " Specified translation : " << separation_
                    << endl;
            }

            transform_ = transformer(separation_);

            break;
        }
        default:
        {
            if (debug)
            {
                Pout<< "patch:" << name()
                    << " Assuming cyclic AMI pairs are colocated" << endl;
            }

            separation_ = Zero;

            transform_ = transformer();

            break;
        }
    }

    if (debug)
    {
        Pout<< "patch: " << name() << nl
            << "    transform = " << transform() << nl << endl;
    }
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::cyclicAMIPolyPatch::resetAMI() const
{
    if (owner())
    {
        const polyPatch& nbr = nbrPatch();
        pointField nbrPoints
        (
            nbrPatch().boundaryMesh().mesh().points(),
            nbrPatch().meshPoints()
        );

        if (debug)
        {
            const Time& t = boundaryMesh().mesh().time();
            OFstream os(t.path()/name() + "_neighbourPatch-org.obj");
            meshTools::writeOBJ(os, nbrPatch().localFaces(), nbrPoints);
        }

        // Transform neighbour patch to local system
        transformPosition(nbrPoints);
        primitivePatch nbrPatch0
        (
            SubList<face>
            (
                nbr.localFaces(),
                nbr.size()
            ),
            nbrPoints
        );

        if (debug)
        {
            const Time& t = boundaryMesh().mesh().time();
            OFstream osN(t.path()/name() + "_neighbourPatch-trans.obj");
            meshTools::writeOBJ(osN, nbrPatch0.localFaces(), nbrPoints);

            OFstream osO(t.path()/name() + "_ownerPatch.obj");
            meshTools::writeOBJ(osO, this->localFaces(), localPoints());
        }

        // Construct/apply AMI interpolation to determine addressing and weights
        AMIs_.resize(1);
        AMIs_.set
        (
            0,
            new AMIInterpolation
            (
                *this,
                nbrPatch0,
                surfPtr(),
                faceAreaIntersect::tmMesh,
                AMIRequireMatch_,
                AMIMethod_,
                AMILowWeightCorrection_,
                AMIReverse_
            )
        );

        AMITransforms_.resize(1, transformer::I);

        if (debug)
        {
            Pout<< "cyclicAMIPolyPatch : " << name()
                << " constructed AMI with " << nl
                << "    " << "srcAddress:" << AMIs_[0].srcAddress().size()
                << nl
                << "    " << "tgAddress :" << AMIs_[0].tgtAddress().size()
                << nl << endl;
        }
    }
}


void Foam::cyclicAMIPolyPatch::calcTransforms()
{
    const cyclicAMIPolyPatch& thisPatch = *this;
    vectorField thisPatchAreas(thisPatch.size());
    forAll(thisPatch, facei)
    {
        thisPatchAreas[facei] = thisPatch[facei].area(thisPatch.points());
    }

    const cyclicAMIPolyPatch& nbrPatch = this->nbrPatch();
    vectorField nbrPatchAreas(nbrPatch.size());
    forAll(nbrPatch, facei)
    {
        nbrPatchAreas[facei] = nbrPatch[facei].area(nbrPatch.points());
    }

    calcTransforms
    (
        thisPatch,
        thisPatch.faceCentres(),
        thisPatchAreas,
        nbrPatch.faceCentres(),
        nbrPatchAreas
    );

    if (debug)
    {
        Pout<< "calcTransforms() : patch: " << name() << nl
            << "    transform = " << transform() << nl << endl;
    }
}


void Foam::cyclicAMIPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    // Clear the invalid AMIs and transforms
    AMIs_.clear();
    AMITransforms_.clear();

    polyPatch::initGeometry(pBufs);
}


void Foam::cyclicAMIPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    calcGeometry
    (
        *this,
        faceCentres(),
        faceAreas(),
        faceCellCentres(),
        nbrPatch().faceCentres(),
        nbrPatch().faceAreas(),
        nbrPatch().faceCellCentres()
    );
}


void Foam::cyclicAMIPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    // Clear the invalid AMIs and transforms
    AMIs_.clear();
    AMITransforms_.clear();

    polyPatch::initMovePoints(pBufs, p);

    // See below. Clear out any local geometry
    primitivePatch::movePoints(p);
}


void Foam::cyclicAMIPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::movePoints(pBufs, p);

    calcTransforms();
}


void Foam::cyclicAMIPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    // Clear the invalid AMIs and transforms
    AMIs_.clear();
    AMITransforms_.clear();

    polyPatch::initUpdateMesh(pBufs);
}


void Foam::cyclicAMIPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    polyPatch::updateMesh(pBufs);
}


void Foam::cyclicAMIPolyPatch::clearGeom()
{
    // Clear the invalid AMIs and transforms
    AMIs_.clear();
    AMITransforms_.clear();

    polyPatch::clearGeom();
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const orderingType ordering,
    const bool AMIRequireMatch,
    const AMIInterpolation::interpolationMethod AMIMethod
)
:
    coupledPolyPatch(name, size, start, index, bm, patchType, ordering),
    nbrPatchName_(word::null),
    nbrPatchID_(-1),
    rotationAxis_(Zero),
    rotationCentre_(point::zero),
    rotationAngleDefined_(false),
    rotationAngle_(0.0),
    AMIs_(),
    AMITransforms_(),
    AMIReverse_(false),
    AMIRequireMatch_(AMIRequireMatch),
    AMILowWeightCorrection_(-1.0),
    AMIMethod_(AMIMethod),
    surfPtr_(nullptr),
    surfDict_(fileName("surface"))
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const bool AMIRequireMatch,
    const AMIInterpolation::interpolationMethod AMIMethod
)
:
    coupledPolyPatch(name, dict, index, bm, patchType),
    cyclicTransform(dict),
    nbrPatchName_(dict.lookupOrDefault<word>("neighbourPatch", "")),
    coupleGroup_(dict),
    nbrPatchID_(-1),
    rotationAxis_(Zero),
    rotationCentre_(point::zero),
    rotationAngleDefined_(false),
    rotationAngle_(0.0),
    AMIs_(),
    AMITransforms_(),
    AMIReverse_(dict.lookupOrDefault<bool>("flipNormals", false)),
    AMIRequireMatch_(AMIRequireMatch),
    AMILowWeightCorrection_(dict.lookupOrDefault("lowWeightCorrection", -1.0)),
    AMIMethod_
    (
        dict.found("method")
      ? AMIInterpolation::wordTointerpolationMethod
        (
            dict.lookup("method")
        )
      : AMIMethod
    ),
    surfPtr_(nullptr),
    surfDict_(dict.subOrEmptyDict("surface"))
{
    if (nbrPatchName_ == word::null && !coupleGroup_.valid())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "No \"neighbourPatch\" or \"coupleGroup\" provided."
            << exit(FatalIOError);
    }

    if (nbrPatchName_ == name)
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Neighbour patch name " << nbrPatchName_
            << " cannot be the same as this patch " << name
            << exit(FatalIOError);
    }

    switch (transformType())
    {
        case ROTATIONAL:
        {
            dict.lookup("rotationAxis") >> rotationAxis_;
            dict.lookup("rotationCentre") >> rotationCentre_;
            if (dict.readIfPresent("rotationAngle", rotationAngle_))
            {
                rotationAngleDefined_ = true;
                rotationAngle_ = degToRad(rotationAngle_);

                if (debug)
                {
                    Info<< "rotationAngle: " << rotationAngle_ << " [rad]"
                        <<  endl;
                }
            }

            scalar magRot = mag(rotationAxis_);
            if (magRot < small)
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "Illegal rotationAxis " << rotationAxis_ << endl
                    << "Please supply a non-zero vector."
                    << exit(FatalIOError);
            }
            rotationAxis_ /= magRot;

            break;
        }
        case TRANSLATIONAL:
        {
            dict.lookup("separation") >> separation_;
            break;
        }
        default:
        {
            // No additional info required
        }
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const cyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    nbrPatchName_(pp.nbrPatchName_),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    rotationAngleDefined_(pp.rotationAngleDefined_),
    rotationAngle_(pp.rotationAngle_),
    AMIs_(),
    AMITransforms_(),
    AMIReverse_(pp.AMIReverse_),
    AMIRequireMatch_(pp.AMIRequireMatch_),
    AMILowWeightCorrection_(pp.AMILowWeightCorrection_),
    AMIMethod_(pp.AMIMethod_),
    surfPtr_(nullptr),
    surfDict_(pp.surfDict_)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const cyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    cyclicTransform(pp),
    nbrPatchName_(nbrPatchName),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    rotationAngleDefined_(pp.rotationAngleDefined_),
    rotationAngle_(pp.rotationAngle_),
    AMIs_(),
    AMITransforms_(),
    AMIReverse_(pp.AMIReverse_),
    AMIRequireMatch_(pp.AMIRequireMatch_),
    AMILowWeightCorrection_(pp.AMILowWeightCorrection_),
    AMIMethod_(pp.AMIMethod_),
    surfPtr_(nullptr),
    surfDict_(pp.surfDict_)
{
    if (nbrPatchName_ == name())
    {
        FatalErrorInFunction
            << "Neighbour patch name " << nbrPatchName_
            << " cannot be the same as this patch " << name()
            << exit(FatalError);
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const cyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, mapAddressing, newStart),
    cyclicTransform(pp),
    nbrPatchName_(pp.nbrPatchName_),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    rotationAngleDefined_(pp.rotationAngleDefined_),
    rotationAngle_(pp.rotationAngle_),
    AMIs_(),
    AMITransforms_(),
    AMIReverse_(pp.AMIReverse_),
    AMIRequireMatch_(pp.AMIRequireMatch_),
    AMILowWeightCorrection_(pp.AMILowWeightCorrection_),
    AMIMethod_(pp.AMIMethod_),
    surfPtr_(nullptr),
    surfDict_(pp.surfDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicAMIPolyPatch::~cyclicAMIPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cyclicAMIPolyPatch::nbrPatchID() const
{
    if (nbrPatchID_ == -1)
    {
        nbrPatchID_ = this->boundaryMesh().findPatchID(nbrPatchName());

        if (nbrPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Illegal neighbourPatch name " << nbrPatchName()
                << nl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        // Check that it is a cyclic AMI patch
        const cyclicAMIPolyPatch& nbrPatch =
            refCast<const cyclicAMIPolyPatch>
            (
                this->boundaryMesh()[nbrPatchID_]
            );

        if (nbrPatch.nbrPatchName() != name())
        {
            WarningInFunction
                << "Patch " << name()
                << " specifies neighbour patch " << nbrPatchName()
                << nl << " but that in return specifies "
                << nbrPatch.nbrPatchName() << endl;
        }
    }

    return nbrPatchID_;
}


bool Foam::cyclicAMIPolyPatch::owner() const
{
    return index() < nbrPatchID();
}


const Foam::cyclicAMIPolyPatch& Foam::cyclicAMIPolyPatch::nbrPatch() const
{
    const polyPatch& pp = this->boundaryMesh()[nbrPatchID()];
    return refCast<const cyclicAMIPolyPatch>(pp);
}


const Foam::autoPtr<Foam::searchableSurface>&
Foam::cyclicAMIPolyPatch::surfPtr() const
{
    const word surfType(surfDict_.lookupOrDefault<word>("type", "none"));

    if (!surfPtr_.valid() && owner() && surfType != "none")
    {
        word surfName(surfDict_.lookupOrDefault("name", name()));

        const polyMesh& mesh = boundaryMesh().mesh();

        surfPtr_ =
            searchableSurface::New
            (
                surfType,
                IOobject
                (
                    surfName,
                    mesh.time().constant(),
                    "triSurface",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                surfDict_
            );
    }

    return surfPtr_;
}


const Foam::PtrList<Foam::AMIInterpolation>&
Foam::cyclicAMIPolyPatch::AMIs() const
{
    if (!owner())
    {
        FatalErrorInFunction
            << "AMI interpolators only available to owner patch"
            << abort(FatalError);
    }

    if (AMIs_.empty())
    {
        resetAMI();
    }

    return AMIs_;
}


const Foam::List<Foam::transformer>&
Foam::cyclicAMIPolyPatch::AMITransforms() const
{
    if (!owner())
    {
        FatalErrorInFunction
            << "AMI transforms only available to owner patch"
            << abort(FatalError);
    }

    if (AMIs_.empty())
    {
        resetAMI();
    }

    return AMITransforms_;
}


bool Foam::cyclicAMIPolyPatch::applyLowWeightCorrection() const
{
    if (owner())
    {
        return AMILowWeightCorrection_ > 0;
    }
    else
    {
        return nbrPatch().AMILowWeightCorrection_ > 0;
    }
}


const Foam::scalarField& Foam::cyclicAMIPolyPatch::weightsSum() const
{
    if (owner())
    {
        return AMIs()[0].srcWeightsSum();
    }
    else
    {
        return nbrPatch().AMIs()[0].tgtWeightsSum();
    }
}


const Foam::scalarField& Foam::cyclicAMIPolyPatch::nbrWeightsSum() const
{
    if (owner())
    {
        return AMIs()[0].tgtWeightsSum();
    }
    else
    {
        return nbrPatch().AMIs()[0].srcWeightsSum();
    }
}


void Foam::cyclicAMIPolyPatch::transformPosition(pointField& l) const
{
    if (transform().rotates())
    {
        if (transformType() == ROTATIONAL)
        {
            l = Foam::transform(transform().R(), l - rotationCentre_)
              + rotationCentre_;
        }
        else
        {
            l = Foam::transform(transform().R(), l);
        }
    }
    else if (transform().translates())
    {
        // transformPosition gets called on the receiving side,
        // separation gets calculated on the sending side so subtract
        l -= transform().t();
    }
}


void Foam::cyclicAMIPolyPatch::transformPosition
(
    point& l,
    const label facei
) const
{
    if (transform().rotates())
    {
        if (transformType() == ROTATIONAL)
        {
            l = Foam::transform(transform().R(), l - rotationCentre_)
              + rotationCentre_;
        }
        else
        {
            l = Foam::transform(transform().R(), l);
        }
    }
    else if (transform().translates())
    {
        l -= transform().t();
    }
}


void Foam::cyclicAMIPolyPatch::transformDirection
(
    vector& d,
    const label facei
) const
{
    if (transform().rotates())
    {
        d = Foam::transform(transform().R(), d);
    }
}


void Foam::cyclicAMIPolyPatch::reverseTransformPosition
(
    point& l,
    const label facei
) const
{
    if (transform().rotates())
    {
        if (transformType() == ROTATIONAL)
        {
            l = Foam::transform(transform().R().T(), l - rotationCentre_)
              + rotationCentre_;
        }
        else
        {
            l = Foam::transform(transform().R().T(), l);
        }
    }
    else if (transform().translates())
    {
        l += transform().t();
    }
}


void Foam::cyclicAMIPolyPatch::reverseTransformDirection
(
    vector& d,
    const label facei
) const
{
    if (transform().rotates())
    {
        d = Foam::transform(transform().R().T(), d);
    }
}


Foam::tmp<Foam::scalarField> Foam::cyclicAMIPolyPatch::interpolate
(
    const scalarField& fld,
    const direction cmpt,
    const direction rank,
    const scalarUList& defaultValues
) const
{
    const cyclicAMIPolyPatch& nei = nbrPatch();

    tmp<scalarField> result(new scalarField(size(), Zero));

    if (owner())
    {
        forAll(AMIs(), i)
        {
            const scalar r =
                pow(inv(AMITransforms()[i]).R()(cmpt, cmpt), rank);

            result.ref() +=
                AMIs()[i].interpolateToSource(r*fld, defaultValues);
        }
    }
    else
    {
        forAll(nei.AMIs(), i)
        {
            const scalar r =
                pow(nei.AMITransforms()[i].R()(cmpt, cmpt), rank);

            result.ref() +=
                nei.AMIs()[i].interpolateToTarget(r*fld, defaultValues);
        }
    }

    return result;
}


void Foam::cyclicAMIPolyPatch::calcGeometry
(
    const primitivePatch& referPatch,
    const pointField& thisCtrs,
    const vectorField& thisAreas,
    const pointField& thisCc,
    const pointField& nbrCtrs,
    const vectorField& nbrAreas,
    const pointField& nbrCc
)
{
    calcTransforms
    (
        referPatch,
        thisCtrs,
        thisAreas,
        nbrCtrs,
        nbrAreas
    );
}


void Foam::cyclicAMIPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{}


bool Foam::cyclicAMIPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    return false;
}


Foam::labelPair Foam::cyclicAMIPolyPatch::pointAMIAndFace
(
    const label facei,
    const vector& n,
    point& p
) const
{
    point pt(p);
    reverseTransformPosition(pt, facei);

    vector nt(n);
    reverseTransformDirection(nt, facei);

    if (owner())
    {
        forAll(AMIs(), i)
        {
            point ptt = AMITransforms()[i].transformPosition(pt);
            const vector ntt = AMITransforms()[i].transform(nt);

            const label nbrFacei =
                AMIs()[i].tgtPointFace(*this, nbrPatch(), ntt, facei, ptt);

            if (nbrFacei >= 0)
            {
                p = ptt;
                return labelPair(i, nbrFacei);
            }
        }
    }
    else
    {
        forAll(nbrPatch().AMIs(), i)
        {
            point ptt =
                nbrPatch().AMITransforms()[i].invTransformPosition(pt);
            const vector ntt =
                nbrPatch().AMITransforms()[i].invTransform(nt);

            const label nbrFacei =
                nbrPatch().AMIs()[i].srcPointFace
                (
                    nbrPatch(),
                    *this,
                    ntt,
                    facei,
                    ptt
                );

            if (nbrFacei >= 0)
            {
                p = ptt;
                return labelPair(i, nbrFacei);
            }
        }
    }

    return labelPair(-1, -1);
}


Foam::label Foam::cyclicAMIPolyPatch::singlePatchProc() const
{
    const cyclicAMIPolyPatch& patch = owner() ? *this : nbrPatch();

    const label proc = patch.AMIs()[0].singlePatchProc();

    for (label i = 1; i < patch.AMIs().size(); ++ i)
    {
        if (patch.AMIs()[i].singlePatchProc() != proc)
        {
            return -1;
        }
    }

    return proc;
}


void Foam::cyclicAMIPolyPatch::write(Ostream& os) const
{
    coupledPolyPatch::write(os);
    if (!nbrPatchName_.empty())
    {
        writeEntry(os, "neighbourPatch", nbrPatchName_);
    }
    coupleGroup_.write(os);

    cyclicTransform::write(os);

    switch (transformType())
    {
        case ROTATIONAL:
        {
            writeEntry(os, "rotationAxis", rotationAxis_);
            writeEntry(os, "rotationCentre", rotationCentre_);

            if (rotationAngleDefined_)
            {
                writeEntry(os, "rotationAngle", radToDeg(rotationAngle_));
            }

            break;
        }
        case TRANSLATIONAL:
        {
            writeEntry(os, "separation", separation_);
            break;
        }
        default:
        {
            // No additional info to write
        }
    }

    if (AMIReverse_)
    {
        writeEntry(os, "flipNormals", AMIReverse_);
    }

    if (AMILowWeightCorrection_ > 0)
    {
        writeEntry(os, "lowWeightCorrection", AMILowWeightCorrection_);
    }

    writeEntry
    (
        os,
        "method",
        AMIInterpolation::interpolationMethodToWord(AMIMethod_)
    );

    if (!surfDict_.empty())
    {
        writeKeyword(os, surfDict_.dictName());
        os  << surfDict_;
    }
}


// ************************************************************************* //
