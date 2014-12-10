/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "streamLineParticle.H"
#include "vectorFieldIOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
//    defineParticleTypeNameAndDebug(streamLineParticle, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::streamLineParticle::calcSubCycleDeltaT
(
    trackingData& td,
    const scalar dt,
    const vector& U
) const
{
    particle testParticle(*this);

    bool oldKeepParticle = td.keepParticle;
    bool oldSwitchProcessor = td.switchProcessor;
    scalar fraction = testParticle.trackToFace(position()+dt*U, td);
    td.keepParticle = oldKeepParticle;
    td.switchProcessor = oldSwitchProcessor;
    // Adapt the dt to subdivide the trajectory into substeps.
    return dt*fraction/td.nSubCycle_;
}


Foam::vector Foam::streamLineParticle::interpolateFields
(
    const trackingData& td,
    const point& position,
    const label cellI,
    const label faceI
)
{
    if (cellI == -1)
    {
        FatalErrorIn("streamLineParticle::interpolateFields(..)")
            << "Cell:" << cellI << abort(FatalError);
    }

    sampledScalars_.setSize(td.vsInterp_.size());
    forAll(td.vsInterp_, scalarI)
    {
        sampledScalars_[scalarI].append
        (
            td.vsInterp_[scalarI].interpolate
            (
                position,
                cellI,
                faceI
            )
        );
    }

    sampledVectors_.setSize(td.vvInterp_.size());
    forAll(td.vvInterp_, vectorI)
    {
        sampledVectors_[vectorI].append
        (
            td.vvInterp_[vectorI].interpolate
            (
                position,
                cellI,
                faceI
            )
        );
    }

    const DynamicList<vector>& U = sampledVectors_[td.UIndex_];

    return U.last();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::streamLineParticle::streamLineParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label lifeTime
)
:
    particle(mesh, position, cellI),
    lifeTime_(lifeTime)
{}


Foam::streamLineParticle::streamLineParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields)
{
    if (readFields)
    {
        //if (is.format() == IOstream::ASCII)
        List<scalarList> sampledScalars;
        List<vectorList> sampledVectors;

        is  >> lifeTime_ >> sampledPositions_ >> sampledScalars
            >> sampledVectors;

        sampledScalars_.setSize(sampledScalars.size());
        forAll(sampledScalars, i)
        {
            sampledScalars_[i].transfer(sampledScalars[i]);
        }
        sampledVectors_.setSize(sampledVectors.size());
        forAll(sampledVectors, i)
        {
            sampledVectors_[i].transfer(sampledVectors[i]);
        }
    }

    // Check state of Istream
    is.check
    (
        "streamLineParticle::streamLineParticle"
        "(const Cloud<streamLineParticle>&, Istream&, bool)"
    );
}


Foam::streamLineParticle::streamLineParticle
(
    const streamLineParticle& p
)
:
    particle(p),
    lifeTime_(p.lifeTime_),
    sampledPositions_(p.sampledPositions_),
    sampledScalars_(p.sampledScalars_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::streamLineParticle::move
(
    trackingData& td,
    const scalar trackTime
)
{
    streamLineParticle& p = static_cast<streamLineParticle&>(*this);

    td.switchProcessor = false;
    td.keepParticle = true;

    scalar tEnd = (1.0 - stepFraction())*trackTime;
    scalar maxDt = mesh_.bounds().mag();

    while
    (
        td.keepParticle
    && !td.switchProcessor
    && lifeTime_ > 0
    )
    {
        // set the lagrangian time-step
        scalar dt = maxDt;

        // Cross cell in steps:
        // - at subiter 0 calculate dt to cross cell in nSubCycle steps
        // - at the last subiter do all of the remaining track
        for (label subIter = 0; subIter < td.nSubCycle_; subIter++)
        {
            --lifeTime_;

            // Store current position and sampled velocity.
            sampledPositions_.append(position());
            vector U = interpolateFields(td, position(), cell(), face());

            if (!td.trackForward_)
            {
                U = -U;
            }

            scalar magU = mag(U);

            if (magU < SMALL)
            {
                // Stagnant particle. Might as well stop
                lifeTime_ = 0;
                break;
            }

            U /= magU;

            if (td.trackLength_ < GREAT)
            {
                dt = td.trackLength_;
                //Pout<< "    subiteration " << subIter
                //    << " : fixed length: updated dt:" << dt << endl;
            }
            else if (subIter == 0 && td.nSubCycle_ > 1)
            {
                // Adapt dt to cross cell in a few steps
                dt = calcSubCycleDeltaT(td, dt, U);
            }
            else if (subIter == td.nSubCycle_ - 1)
            {
                // Do full step on last subcycle
                dt = maxDt;
            }


            scalar fraction = trackToFace(position() + dt*U, td);
            dt *= fraction;

            tEnd -= dt;
            stepFraction() = 1.0 - tEnd/trackTime;

            if (tEnd <= ROOTVSMALL)
            {
                // Force removal
                lifeTime_ = 0;
            }

            if
            (
                face() != -1
            || !td.keepParticle
            ||  td.switchProcessor
            ||  lifeTime_ == 0
            )
            {
                break;
            }
        }
    }


    if (!td.keepParticle || lifeTime_ == 0)
    {
        if (lifeTime_ == 0)
        {
            if (debug)
            {
                Pout<< "streamLineParticle : Removing stagnant particle:"
                    << p.position()
                    << " sampled positions:" << sampledPositions_.size()
                    << endl;
            }
            td.keepParticle = false;
        }
        else
        {
            // Normal exit. Store last position and fields
            sampledPositions_.append(position());
            interpolateFields(td, position(), cell(), face());

            if (debug)
            {
                Pout<< "streamLineParticle : Removing particle:"
                    << p.position()
                    << " sampled positions:" << sampledPositions_.size()
                    << endl;
            }
        }

        // Transfer particle data into trackingData.
        //td.allPositions_.append(sampledPositions_);
        td.allPositions_.append(vectorList());
        vectorList& top = td.allPositions_.last();
        top.transfer(sampledPositions_);

        forAll(sampledScalars_, i)
        {
            //td.allScalars_[i].append(sampledScalars_[i]);
            td.allScalars_[i].append(scalarList());
            scalarList& top = td.allScalars_[i].last();
            top.transfer(sampledScalars_[i]);
        }
        forAll(sampledVectors_, i)
        {
            //td.allVectors_[i].append(sampledVectors_[i]);
            td.allVectors_[i].append(vectorList());
            vectorList& top = td.allVectors_[i].last();
            top.transfer(sampledVectors_[i]);
        }
    }

    return td.keepParticle;
}


bool Foam::streamLineParticle::hitPatch
(
    const polyPatch&,
    trackingData& td,
    const label patchI,
    const scalar trackFraction,
    const tetIndices& tetIs
)
{
    // Disable generic patch interaction
    return false;
}


void Foam::streamLineParticle::hitWedgePatch
(
    const wedgePolyPatch& pp,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitSymmetryPlanePatch
(
    const symmetryPlanePolyPatch& pp,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitSymmetryPatch
(
    const symmetryPolyPatch& pp,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitCyclicPatch
(
    const cyclicPolyPatch& pp,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitProcessorPatch
(
    const processorPolyPatch&,
    trackingData& td
)
{
    // Switch particle
    td.switchProcessor = true;
}


void Foam::streamLineParticle::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackingData& td,
    const tetIndices&
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitPatch
(
    const polyPatch& wpp,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::readFields(Cloud<streamLineParticle>& c)
{
    if (!c.size())
    {
        return;
    }

    particle::readFields(c);

    IOField<label> lifeTime
    (
        c.fieldIOobject("lifeTime", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, lifeTime);

    vectorFieldIOField sampledPositions
    (
        c.fieldIOobject("sampledPositions", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, sampledPositions);

//    vectorFieldIOField sampleVelocity
//    (
//        c.fieldIOobject("sampleVelocity", IOobject::MUST_READ)
//    );
//    c.checkFieldIOobject(c, sampleVelocity);

    label i = 0;
    forAllIter(Cloud<streamLineParticle>, c, iter)
    {
        iter().lifeTime_ = lifeTime[i];
        iter().sampledPositions_.transfer(sampledPositions[i]);
//        iter().sampleVelocity_.transfer(sampleVelocity[i]);
        i++;
    }
}


void Foam::streamLineParticle::writeFields(const Cloud<streamLineParticle>& c)
{
    particle::writeFields(c);

    label np =  c.size();

    IOField<label> lifeTime
    (
        c.fieldIOobject("lifeTime", IOobject::NO_READ),
        np
    );
    vectorFieldIOField sampledPositions
    (
        c.fieldIOobject("sampledPositions", IOobject::NO_READ),
        np
    );
//    vectorFieldIOField sampleVelocity
//    (
//        c.fieldIOobject("sampleVelocity", IOobject::NO_READ),
//        np
//    );

    label i = 0;
    forAllConstIter(Cloud<streamLineParticle>, c, iter)
    {
        lifeTime[i] = iter().lifeTime_;
        sampledPositions[i] = iter().sampledPositions_;
//        sampleVelocity[i] = iter().sampleVelocity_;
        i++;
    }

    lifeTime.write();
    sampledPositions.write();
//    sampleVelocity.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const streamLineParticle& p)
{
    os  << static_cast<const particle&>(p)
        << token::SPACE << p.lifeTime_
        << token::SPACE << p.sampledPositions_
        << token::SPACE << p.sampledScalars_
        << token::SPACE << p.sampledVectors_;

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const streamLineParticle&)");

    return os;
}


// ************************************************************************* //
