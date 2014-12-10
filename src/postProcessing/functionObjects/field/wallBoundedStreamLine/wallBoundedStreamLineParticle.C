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

#include "wallBoundedStreamLineParticle.H"
#include "vectorFieldIOField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector Foam::wallBoundedStreamLineParticle::interpolateFields
(
    const trackingData& td,
    const point& position,
    const label cellI,
    const label faceI
)
{
    if (cellI == -1)
    {
        FatalErrorIn("wallBoundedStreamLineParticle::interpolateFields(..)")
            << "Cell:" << cellI << abort(FatalError);
    }

    const tetIndices ti = currentTetIndices();

    const vector U = td.vvInterp_[td.UIndex_].interpolate
    (
        position,
        ti,     //cellI,
        faceI
    );

    // Check if at different position
    if
    (
       !sampledPositions_.size()
     || magSqr(sampledPositions_.last()-position) > Foam::sqr(SMALL)
    )
    {
        // Store the location
        sampledPositions_.append(position);

        // Store the scalar fields
        sampledScalars_.setSize(td.vsInterp_.size());
        forAll(td.vsInterp_, scalarI)
        {
            sampledScalars_[scalarI].append
            (
                td.vsInterp_[scalarI].interpolate
                (
                    position,
                    ti,     //cellI,
                    faceI
                )
            );
        }

        // Store the vector fields
        sampledVectors_.setSize(td.vvInterp_.size());
        forAll(td.vvInterp_, vectorI)
        {
            vector positionU;
            if (vectorI == td.UIndex_)
            {
                positionU = U;
            }
            else
            {
                positionU = td.vvInterp_[vectorI].interpolate
                (
                    position,
                    ti,     //cellI,
                    faceI
                );
            }
            sampledVectors_[vectorI].append(positionU);
        }
    }

    return U;
}


Foam::vector Foam::wallBoundedStreamLineParticle::sample
(
    trackingData& td
)
{
    vector U = interpolateFields(td, position(), cell(), tetFace());

    if (!td.trackForward_)
    {
        U = -U;
    }

    scalar magU = mag(U);

    if (magU < SMALL)
    {
        // Stagnant particle. Might as well stop
        lifeTime_ = 0;
    }

    U /= magU;

    return U;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoundedStreamLineParticle::wallBoundedStreamLineParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label meshEdgeStart,
    const label diagEdge,
    const label lifeTime
)
:
    wallBoundedParticle
    (
        mesh,
        position,
        cellI,
        tetFaceI,
        tetPtI,
        meshEdgeStart,
        diagEdge
    ),
    lifeTime_(lifeTime)
{}


Foam::wallBoundedStreamLineParticle::wallBoundedStreamLineParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    wallBoundedParticle(mesh, is, readFields)
{
    if (readFields)
    {
        List<scalarList> sampledScalars;
        List<vectorList> sampledVectors;

        is  >> lifeTime_
            >> sampledPositions_ >> sampledScalars >> sampledVectors;

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
        "wallBoundedStreamLineParticle::wallBoundedStreamLineParticle"
        "(const Cloud<wallBoundedStreamLineParticle>&, Istream&, bool)"
    );
}


Foam::wallBoundedStreamLineParticle::wallBoundedStreamLineParticle
(
    const wallBoundedStreamLineParticle& p
)
:
    wallBoundedParticle(p),
    lifeTime_(p.lifeTime_),
    sampledPositions_(p.sampledPositions_),
    sampledScalars_(p.sampledScalars_),
    sampledVectors_(p.sampledVectors_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::wallBoundedStreamLineParticle::move
(
    trackingData& td,
    const scalar trackTime
)
{
    wallBoundedStreamLineParticle& p = static_cast
    <
        wallBoundedStreamLineParticle&
    >(*this);


    // Check position is inside tet
    //checkInside();

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

        --lifeTime_;

        // Get sampled velocity and fields. Store if position changed.
        vector U = sample(td);

        // !user parameter!
        if (dt < SMALL)
        {
            // Force removal
            lifeTime_ = 0;
            break;
        }


        if (td.trackLength_ < GREAT)
        {
            dt = td.trackLength_;
        }


        scalar fraction = trackToEdge(td, position() + dt*U);
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
            !td.keepParticle
        ||  td.switchProcessor
        ||  lifeTime_ == 0
        )
        {
            break;
        }
    }


    if (!td.keepParticle || lifeTime_ == 0)
    {
        if (lifeTime_ == 0)
        {
            if (debug)
            {
                Pout<< "wallBoundedStreamLineParticle :"
                    << " Removing stagnant particle:"
                    << p.position()
                    << " sampled positions:" << sampledPositions_.size()
                    << endl;
            }
            td.keepParticle = false;
        }
        else
        {
            // Normal exit. Store last position and fields
            sample(td);

            if (debug)
            {
                Pout<< "wallBoundedStreamLineParticle : Removing particle:"
                    << p.position()
                    << " sampled positions:" << sampledPositions_.size()
                    << endl;
            }
        }

        // Transfer particle data into trackingData.
        {
            //td.allPositions_.append(sampledPositions_);
            td.allPositions_.append(vectorList());
            vectorList& top = td.allPositions_.last();
            top.transfer(sampledPositions_);
        }

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


void Foam::wallBoundedStreamLineParticle::readFields
(
    Cloud<wallBoundedStreamLineParticle>& c
)
{
    if (!c.size())
    {
        return;
    }

    wallBoundedParticle::readFields(c);

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

    label i = 0;
    forAllIter(Cloud<wallBoundedStreamLineParticle>, c, iter)
    {
        iter().lifeTime_ = lifeTime[i];
        iter().sampledPositions_.transfer(sampledPositions[i]);
        i++;
    }
}


void Foam::wallBoundedStreamLineParticle::writeFields
(
    const Cloud<wallBoundedStreamLineParticle>& c
)
{
    wallBoundedParticle::writeFields(c);

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

    label i = 0;
    forAllConstIter(Cloud<wallBoundedStreamLineParticle>, c, iter)
    {
        lifeTime[i] = iter().lifeTime_;
        sampledPositions[i] = iter().sampledPositions_;
        i++;
    }

    lifeTime.write();
    sampledPositions.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const wallBoundedStreamLineParticle& p
)
{
    os  << static_cast<const wallBoundedParticle&>(p)
        << token::SPACE << p.lifeTime_
        << token::SPACE << p.sampledPositions_
        << token::SPACE << p.sampledScalars_
        << token::SPACE << p.sampledVectors_;

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const wallBoundedStreamLineParticle&)"
    );

    return os;
}


// ************************************************************************* //
