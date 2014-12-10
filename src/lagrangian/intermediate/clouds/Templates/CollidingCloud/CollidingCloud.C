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

#include "CollidingCloud.H"
#include "CollisionModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::CollidingCloud<CloudType>::setModels()
{
    collisionModel_.reset
    (
        CollisionModel<CollidingCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}


template<class CloudType>
template<class TrackData>
void  Foam::CollidingCloud<CloudType>::moveCollide
(
    TrackData& td,
    const scalar deltaT
)
{
    td.part() = TrackData::tpVelocityHalfStep;
    CloudType::move(td,  deltaT);

    td.part() = TrackData::tpLinearTrack;
    CloudType::move(td,  deltaT);

    // td.part() = TrackData::tpRotationalTrack;
    // CloudType::move(td);

    this->updateCellOccupancy();

    this->collision().collide();

    td.part() = TrackData::tpVelocityHalfStep;
    CloudType::move(td,  deltaT);
}



template<class CloudType>
void Foam::CollidingCloud<CloudType>::cloudReset(CollidingCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    collisionModel_.reset(c.collisionModel_.ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CollidingCloud<CloudType>::CollidingCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g,
    bool readFields
)
:
    CloudType(cloudName, rho, U, mu, g, false),
    constProps_(this->particleProperties()),
    collisionModel_(NULL)
{
    if (this->solution().steadyState())
    {
        FatalErrorIn
        (
            "Foam::CollidingCloud<CloudType>::CollidingCloud"
            "("
                "const word&, "
                "const volScalarField&, "
                "const volVectorField&, "
                "const volScalarField&, "
                "const dimensionedVector&, "
                "bool"
            ")"
        )   << "Collision modelling not currently available for steady state "
            << "calculations" << exit(FatalError);
    }

    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this);
        }
    }
}


template<class CloudType>
Foam::CollidingCloud<CloudType>::CollidingCloud
(
    CollidingCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    collisionModel_(c.collisionModel_->clone())
{}


template<class CloudType>
Foam::CollidingCloud<CloudType>::CollidingCloud
(
    const fvMesh& mesh,
    const word& name,
    const CollidingCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    collisionModel_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CollidingCloud<CloudType>::~CollidingCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::CollidingCloud<CloudType>::hasWallImpactDistance() const
{
    return !collision().controlsWallInteraction();
}


template<class CloudType>
void Foam::CollidingCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<CollidingCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::CollidingCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::CollidingCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::template
            TrackingData<CollidingCloud<CloudType> > td(*this);

        this->solve(td);
    }
}


template<class CloudType>
template<class TrackData>
void  Foam::CollidingCloud<CloudType>::motion(TrackData& td)
{
    // Sympletic leapfrog integration of particle forces:
    // + apply half deltaV with stored force
    // + move positions with new velocity
    // + calculate forces in new position
    // + apply half deltaV with new force

    label nSubCycles = collision().nSubCycles();

    if (nSubCycles > 1)
    {
        Info<< "    " << nSubCycles << " move-collide subCycles" << endl;

        subCycleTime moveCollideSubCycle
        (
            const_cast<Time&>(this->db().time()),
            nSubCycles
        );

        while(!(++moveCollideSubCycle).end())
        {
            moveCollide(td, this->db().time().deltaTValue());
        }

        moveCollideSubCycle.endSubCycle();
    }
    else
    {
        moveCollide(td, this->db().time().deltaTValue());
    }
}


template<class CloudType>
void Foam::CollidingCloud<CloudType>::info()
{
    CloudType::info();

    scalar rotationalKineticEnergy = rotationalKineticEnergyOfSystem();
    reduce(rotationalKineticEnergy, sumOp<scalar>());

    Info<< "    Rotational kinetic energy       = "
        << rotationalKineticEnergy << nl;
}


// ************************************************************************* //
