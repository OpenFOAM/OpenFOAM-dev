/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "streamlinesParticle.H"
#include "streamlinesCloud.H"
#include "vectorFieldIOField.H"
#include "scalarFieldIOField.H"
#include "transformerIOList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector Foam::streamlinesParticle::interpolateFields
(
    const trackingData& td,
    const point& position,
    const label celli,
    const label facei
)
{
    if (celli == -1)
    {
        FatalErrorInFunction
            << "Cell:" << celli << abort(FatalError);
    }

    bool interpolatedU = false;
    vector U = vector::uniform(NaN);

    forAll(td.scalarInterp_, fieldi)
    {
        #define InterpolateType(Type, nullArg)                        \
            if (td.Type##Interp_.set(fieldi))                         \
            {                                                         \
                const Type s =                                        \
                    td.Type##Interp_[fieldi].interpolate              \
                    (                                                 \
                        position,                                     \
                        celli,                                        \
                        facei                                         \
                    );                                                \
                                                                      \
                sampled##Type##s_.setSize(td.Type##Interp_.size());   \
                sampled##Type##s_[fieldi].append                      \
                (                                                     \
                    td.trackOutside_ ? transform_.invTransform(s) : s \
                );                                                    \
            }
        FOR_ALL_FIELD_TYPES(InterpolateType);
        #undef InterpolateType

        if
        (
            td.vectorInterp_.set(fieldi)
         && &td.vectorInterp_[fieldi] == &td.UInterp_
        )
        {
            interpolatedU = true;
            U = sampledvectors_[fieldi].last();
        }
    }

    // Interpolate the velocity if it has not already been done
    if (!interpolatedU)
    {
        U = td.UInterp_.interpolate(position, celli, facei);
    }

    return U;
}


void Foam::streamlinesParticle::endTrack(trackingData& td)
{
    const label n = sampledPositions_.size();

    const label trackPartIndex =
        td.trackForward_ ? trackPartIndex_ : -1 - trackPartIndex_;

    if (!td.trackForward_) reverse(sampledPositions_);
    td.allPositions_.append(sampledPositions_);
    sampledPositions_.clear();

    td.allTracks_.append(List<label>(n, trackIndex_));
    td.allTrackParts_.append(List<label>(n, trackPartIndex));
    trackPartIndex_ ++;

    if (!td.trackForward_) reverse(sampledAges_);
    td.allAges_.append(sampledAges_);
    sampledAges_.clear();

    forAll(td.scalarInterp_, fieldi)
    {
        #define EndTrackType(Type, nullArg)                                 \
            if (td.Type##Interp_.set(fieldi))                               \
            {                                                               \
                if (!td.trackForward_) reverse(sampled##Type##s_[fieldi]);  \
                td.all##Type##s_[fieldi].append(sampled##Type##s_[fieldi]); \
                sampled##Type##s_[fieldi].clear();                          \
            }
        FOR_ALL_FIELD_TYPES(EndTrackType);
        #undef EndTrackType
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::streamlinesParticle::streamlinesParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    label& nLocateBoundaryHits,
    const label lifeTime,
    const label trackIndex
)
:
    particle(mesh, position, celli, nLocateBoundaryHits),
    lifeTime_(lifeTime),
    trackIndex_(trackIndex),
    trackPartIndex_(0),
    age_(0),
    transform_(transformer::I)
{}


Foam::streamlinesParticle::streamlinesParticle(Istream& is, bool readFields)
:
    particle(is, readFields)
{
    if (readFields)
    {
        is  >> lifeTime_ >> trackIndex_ >> trackPartIndex_ >> age_
            >> transform_ >> sampledPositions_ >> sampledAges_;

        #define ReadSampledTypes(Type, nullArg)                     \
            List<List<Type>> sampled##Type##s;                      \
            is >> sampled##Type##s;                                 \
            sampled##Type##s_.setSize(sampled##Type##s.size());     \
            forAll(sampled##Type##s, i)                             \
            {                                                       \
                sampled##Type##s_[i].transfer(sampled##Type##s[i]); \
            }
        FOR_ALL_FIELD_TYPES(ReadSampledTypes);
        #undef ReadSampledTypes
    }

    // Check state of Istream
    is.check
    (
        "streamlinesParticle::streamlinesParticle"
        "(const Cloud<streamlinesParticle>&, Istream&, bool)"
    );
}


Foam::streamlinesParticle::streamlinesParticle
(
    const streamlinesParticle& p
)
:
    particle(p),
    lifeTime_(p.lifeTime_),
    trackIndex_(p.trackIndex_),
    trackPartIndex_(p.trackPartIndex_),
    age_(p.age_),
    transform_(p.transform_),
    sampledPositions_(p.sampledPositions_),
    sampledAges_(p.sampledAges_)
    #define SampledTypesInit(Type, nullArg) \
        , sampled##Type##s_(p.sampled##Type##s_)
    FOR_ALL_FIELD_TYPES(SampledTypesInit)
    #undef SampledTypesInit
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::streamlinesParticle::move
(
    streamlinesCloud& cloud,
    trackingData& td
)
{
    td.keepParticle = true;
    td.sendToProc = -1;

    const scalar maxDt = td.mesh.bounds().mag();

    while (td.keepParticle && td.sendToProc == -1 && lifeTime_ > 0)
    {
        scalar dt = maxDt;

        // Cross cell in steps:
        // - at subiter 0 calculate dt to cross cell in nSubCycle steps
        // - at the last subiter do all of the remaining track
        for (label subIter = 0; subIter < max(1, td.nSubCycle_); subIter++)
        {
            --lifeTime_;

            // Store current position and sampled velocity.
            sampledPositions_.append
            (
                td.trackOutside_
              ? transform_.invTransformPosition(position(td.mesh))
              : position(td.mesh)
            );
            sampledAges_.append(age_);
            vector U = interpolateFields(td, position(td.mesh), cell(), face());

            if (!td.trackForward_)
            {
                U = -U;
            }

            scalar magU = mag(U);

            if (magU < small)
            {
                // Stagnant particle. Might as well stop
                lifeTime_ = 0;
                break;
            }

            U /= magU;

            if (td.trackLength_ < great)
            {
                // No sub-cycling. Track a set length on each step.
                dt = td.trackLength_;
            }
            else if (subIter == 0)
            {
                // Sub-cycling. Cross the cell in nSubCycle steps.
                particle copy(*this);
                copy.trackToFace(td.mesh, maxDt*U, 1);
                dt *= (copy.stepFraction() - stepFraction())/td.nSubCycle_;
            }
            else if (subIter == td.nSubCycle_ - 1)
            {
                // Sub-cycling. Track the whole cell on the last step.
                dt = maxDt;
            }

            age_ +=
                (td.trackForward_ ? +1 : -1)
               *dt/magU
               *(1 - trackToAndHitFace(dt*U, 0, cloud, td));

            if (!td.keepParticle || td.sendToProc != -1 || lifeTime_ == 0)
            {
                break;
            }
        }
    }

    if (!td.keepParticle || lifeTime_ == 0)
    {
        if (lifeTime_ == 0)
        {
            // Failure exit. Particle stagnated or its life ran out.
            if (debug)
            {
                Pout<< "streamlinesParticle: Removing stagnant particle:"
                    << position(td.mesh) << " sampled positions:"
                    << sampledPositions_.size() << endl;
            }
            td.keepParticle = false;
        }
        else
        {
            // Normal exit. Store last position and fields
            sampledPositions_.append
            (
                td.trackOutside_
              ? transform_.invTransformPosition(position(td.mesh))
              : position(td.mesh)
            );
            sampledAges_.append(age_);
            interpolateFields(td, position(td.mesh), cell(), face());

            if (debug)
            {
                Pout<< "streamlinesParticle: Removing particle:"
                    << position(td.mesh) << " sampled positions:"
                    << sampledPositions_.size() << endl;
            }
        }

        // End this track
        endTrack(td);
    }

    return td.keepParticle;
}


void Foam::streamlinesParticle::hitWedgePatch
(
    streamlinesCloud&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamlinesParticle::hitSymmetryPlanePatch
(
    streamlinesCloud&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamlinesParticle::hitSymmetryPatch
(
    streamlinesCloud&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamlinesParticle::hitCyclicPatch
(
    streamlinesCloud& cloud,
    trackingData& td
)
{
    const cyclicPolyPatch& cpp =
        static_cast<const cyclicPolyPatch&>
        (
            td.mesh.boundaryMesh()[patch(td.mesh)]
        );

    // End this track
    if (!td.trackOutside_ && cpp.transform().transformsPosition())
    {
        endTrack(td);
    }

    // Move across the cyclic
    particle::hitCyclicPatch(cloud, td);
}


bool Foam::streamlinesParticle::hitNonConformalCyclicPatch
(
    const vector& displacement,
    const scalar fraction,
    const label patchi,
    streamlinesCloud& cloud,
    trackingData& td
)
{
    // Move across the cyclic
    const bool result =
        particle::hitNonConformalCyclicPatch
        (
            displacement,
            fraction,
            patchi,
            cloud,
            td
        );

    // End this track
    if (result) endTrack(td);

    return result;
}

void Foam::streamlinesParticle::hitProcessorPatch
(
    streamlinesCloud& cloud,
    trackingData& td
)
{
    const processorPolyPatch& ppp =
        static_cast<const processorPolyPatch&>
        (
            td.mesh.boundaryMesh()[patch(td.mesh)]
        );

    // End this track
    if (!td.trackOutside_ && ppp.transform().transformsPosition())
    {
        endTrack(td);
    }

    particle::hitProcessorPatch(cloud, td);
}


void Foam::streamlinesParticle::hitWallPatch
(
    streamlinesCloud&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamlinesParticle::transformProperties
(
    const transformer& transform
)
{
    transform_ = transform & transform_;
}


void Foam::streamlinesParticle::readFields(Cloud<streamlinesParticle>& c)
{
    bool valid = c.size();

    particle::readFields(c);

    IOField<label> lifeTime
    (
        c.fieldIOobject("lifeTime", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, lifeTime);

    IOField<label> trackIndex
    (
        c.fieldIOobject("trackIndex", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, trackIndex);

    IOField<label> trackPartIndex
    (
        c.fieldIOobject("trackPartIndex", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, trackPartIndex);

    IOField<scalar> age
    (
        c.fieldIOobject("age", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, age);

    transformerIOList transform
    (
        c.fieldIOobject("transform", IOobject::MUST_READ),
        valid
    );
    //c.checkFieldIOobject(c, transform);

    vectorFieldIOField sampledPositions
    (
        c.fieldIOobject("sampledPositions", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, sampledPositions);

    scalarFieldIOField sampledAges
    (
        c.fieldIOobject("sampledAges", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, sampledAges);

    label i = 0;
    forAllIter(Cloud<streamlinesParticle>, c, iter)
    {
        iter().lifeTime_ = lifeTime[i];
        iter().trackIndex_ = trackIndex[i];
        iter().trackPartIndex_ = trackPartIndex[i];
        iter().age_ = age[i];
        iter().transform_ = transform[i];
        iter().sampledPositions_.transfer(sampledPositions[i]);
        iter().sampledAges_.transfer(sampledAges[i]);
        i++;
    }
}


void Foam::streamlinesParticle::writeFields(const Cloud<streamlinesParticle>& c)
{
    particle::writeFields(c);

    label np = c.size();

    IOField<label> lifeTime
    (
        c.fieldIOobject("lifeTime", IOobject::NO_READ),
        np
    );

    IOList<label> trackIndex
    (
        c.fieldIOobject("trackIndex", IOobject::NO_READ),
        np
    );

    IOList<label> trackPartIndex
    (
        c.fieldIOobject("trackPartIndex", IOobject::NO_READ),
        np
    );

    IOField<scalar> age
    (
        c.fieldIOobject("age", IOobject::NO_READ),
        np
    );

    transformerIOList transform
    (
        c.fieldIOobject("transform", IOobject::NO_READ),
        np
    );

    vectorFieldIOField sampledPositions
    (
        c.fieldIOobject("sampledPositions", IOobject::NO_READ),
        np
    );

    scalarFieldIOField sampledAges
    (
        c.fieldIOobject("sampledAges", IOobject::NO_READ),
        np
    );

    label i = 0;
    forAllConstIter(Cloud<streamlinesParticle>, c, iter)
    {
        lifeTime[i] = iter().lifeTime_;
        trackIndex[i] = iter().trackIndex_;
        trackPartIndex[i] = iter().trackPartIndex_;
        age[i] = iter().age_;
        transform[i] = iter().transform_;
        sampledPositions[i] = iter().sampledPositions_;
        sampledAges[i] = iter().sampledAges_;
        i++;
    }

    lifeTime.write(np > 0);
    trackIndex.write(np > 0);
    trackPartIndex.write(np > 0);
    age.write(np > 0);
    transform.write(np > 0);
    sampledPositions.write(np > 0);
    sampledAges.write(np > 0);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const streamlinesParticle& p)
{

    os  << static_cast<const particle&>(p)
        << token::SPACE << p.lifeTime_
        << token::SPACE << p.trackIndex_
        << token::SPACE << p.trackPartIndex_
        << token::SPACE << p.age_
        << token::SPACE << p.transform_
        << token::SPACE << p.sampledPositions_
        << token::SPACE << p.sampledAges_
        #define WriteSampledTypes(Type, nullArg) \
            << token::SPACE << p.sampled##Type##s_
        FOR_ALL_FIELD_TYPES(WriteSampledTypes);
        #undef WriteSampledTypes

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const streamlinesParticle&)");

    return os;
}


// ************************************************************************* //
