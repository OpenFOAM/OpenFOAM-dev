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

#include "streamLineParticle.H"
#include "streamLineParticleCloud.H"
#include "vectorFieldIOField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector Foam::streamLineParticle::interpolateFields
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

    sampledScalars_.setSize(td.vsInterp_.size());
    forAll(td.vsInterp_, scalarI)
    {
        sampledScalars_[scalarI].append
        (
            td.vsInterp_[scalarI].interpolate
            (
                position,
                celli,
                facei
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
                celli,
                facei
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
    const label celli,
    const label lifeTime
)
:
    particle(mesh, position, celli),
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
    streamLineParticleCloud& cloud,
    trackingData& td,
    const scalar
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    const scalar maxDt = mesh().bounds().mag();

    while (td.keepParticle && !td.switchProcessor && lifeTime_ > 0)
    {
        scalar dt = maxDt;

        // Cross cell in steps:
        // - at subiter 0 calculate dt to cross cell in nSubCycle steps
        // - at the last subiter do all of the remaining track
        for (label subIter = 0; subIter < max(1, td.nSubCycle_); subIter++)
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
                copy.trackToFace(maxDt*U, 1);
                dt *= (copy.stepFraction() - stepFraction())/td.nSubCycle_;
            }
            else if (subIter == td.nSubCycle_ - 1)
            {
                // Sub-cycling. Track the whole cell on the last step.
                dt = maxDt;
            }

            trackToAndHitFace(dt*U, 0, cloud, td);

            if
            (
                onFace()
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
            // Failure exit. Particle stagnated or it's life ran out.
            if (debug)
            {
                Pout<< "streamLineParticle: Removing stagnant particle:"
                    << position() << " sampled positions:"
                    << sampledPositions_.size() << endl;
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
                Pout<< "streamLineParticle: Removing particle:" << position()
                    << " sampled positions:" << sampledPositions_.size()
                    << endl;
            }
        }

        // Transfer particle data into trackingData.
        td.allPositions_.append(vectorList());
        vectorList& top = td.allPositions_.last();
        top.transfer(sampledPositions_);

        forAll(sampledScalars_, i)
        {
            td.allScalars_[i].append(scalarList());
            scalarList& top = td.allScalars_[i].last();
            top.transfer(sampledScalars_[i]);
        }
        forAll(sampledVectors_, i)
        {
            td.allVectors_[i].append(vectorList());
            vectorList& top = td.allVectors_[i].last();
            top.transfer(sampledVectors_[i]);
        }
    }

    return td.keepParticle;
}


bool Foam::streamLineParticle::hitPatch(streamLineParticleCloud&, trackingData&)
{
    // Disable generic patch interaction
    return false;
}


void Foam::streamLineParticle::hitWedgePatch
(
    streamLineParticleCloud&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitSymmetryPlanePatch
(
    streamLineParticleCloud&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitSymmetryPatch
(
    streamLineParticleCloud&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitCyclicPatch
(
    streamLineParticleCloud&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitCyclicAMIPatch
(
    const vector&,
    const scalar,
    streamLineParticleCloud&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitCyclicACMIPatch
(
    const vector&,
    const scalar,
    streamLineParticleCloud&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitCyclicRepeatAMIPatch
(
    const vector&,
    const scalar,
    streamLineParticleCloud&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitProcessorPatch
(
    streamLineParticleCloud&,
    trackingData& td
)
{
    // Switch particle
    td.switchProcessor = true;
}


void Foam::streamLineParticle::hitWallPatch
(
    streamLineParticleCloud&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::readFields(Cloud<streamLineParticle>& c)
{
//    if (!c.size())
//    {
//        return;
//    }
    bool valid = c.size();

    particle::readFields(c);

    IOField<label> lifeTime
    (
        c.fieldIOobject("lifeTime", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, lifeTime);

    vectorFieldIOField sampledPositions
    (
        c.fieldIOobject("sampledPositions", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, sampledPositions);

    label i = 0;
    forAllIter(Cloud<streamLineParticle>, c, iter)
    {
        iter().lifeTime_ = lifeTime[i];
        iter().sampledPositions_.transfer(sampledPositions[i]);
        i++;
    }
}


void Foam::streamLineParticle::writeFields(const Cloud<streamLineParticle>& c)
{
    particle::writeFields(c);

    label np = c.size();

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

    label i = 0;
    forAllConstIter(Cloud<streamLineParticle>, c, iter)
    {
        lifeTime[i] = iter().lifeTime_;
        sampledPositions[i] = iter().sampledPositions_;
        i++;
    }

    lifeTime.write(np > 0);
    sampledPositions.write(np > 0);
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
