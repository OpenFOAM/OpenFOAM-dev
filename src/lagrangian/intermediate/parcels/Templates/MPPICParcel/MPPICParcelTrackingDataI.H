/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "AveragingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ParcelType>
template<class CloudType>
inline Foam::MPPICParcel<ParcelType>::TrackingData<CloudType>::TrackingData
(
    CloudType& cloud,
    trackPart part
)
:
    ParcelType::template TrackingData<CloudType>(cloud),
    volumeAverage_
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":volumeAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.solution().dict(),
            cloud.mesh()
        )
    ),
    radiusAverage_
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":radiusAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.solution().dict(),
            cloud.mesh()
        )
    ),
    rhoAverage_
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":rhoAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.solution().dict(),
            cloud.mesh()
        )
    ),
    uAverage_
    (
        AveragingMethod<vector>::New
        (
            IOobject
            (
                cloud.name() + ":uAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.solution().dict(),
            cloud.mesh()
        )
    ),
    uSqrAverage_
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":uSqrAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.solution().dict(),
            cloud.mesh()
        )
    ),
    frequencyAverage_
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":frequencyAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.solution().dict(),
            cloud.mesh()
        )
    ),
    massAverage_
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":massAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.solution().dict(),
            cloud.mesh()
        )
    ),
    part_(part)
{}


template<class ParcelType>
template<class CloudType>
inline void
Foam::MPPICParcel<ParcelType>::TrackingData<CloudType>::updateAverages
(
    CloudType& cloud
)
{
    // zero the sums
    volumeAverage_() = 0;
    radiusAverage_() = 0;
    rhoAverage_() = 0;
    uAverage_() = vector::zero;
    uSqrAverage_() = 0;
    frequencyAverage_() = 0;
    massAverage_() = 0;

    // temporary weights
    autoPtr<AveragingMethod<scalar> > weightAveragePtr
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":weightAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.solution().dict(),
            cloud.mesh()
        )
    );
    AveragingMethod<scalar>& weightAverage = weightAveragePtr();

    // averaging sums
    forAllConstIter(typename CloudType, cloud, iter)
    {
        const typename CloudType::parcelType& p = iter();
        const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), cloud.mesh());

        const scalar m = p.nParticle()*p.mass();

        volumeAverage_->add(p.position(), tetIs, p.nParticle()*p.volume());
        rhoAverage_->add(p.position(), tetIs, m*p.rho());
        uAverage_->add(p.position(), tetIs, m*p.U());
        massAverage_->add(p.position(), tetIs, m);
    }
    volumeAverage_->average();
    massAverage_->average();
    rhoAverage_->average(massAverage_);
    uAverage_->average(massAverage_);

    // squared velocity deviation
    forAllConstIter(typename CloudType, cloud, iter)
    {
        const typename CloudType::parcelType& p = iter();
        const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), cloud.mesh());

        const vector u = uAverage_->interpolate(p.position(), tetIs);

        uSqrAverage_->add
        (
            p.position(),
            tetIs,
            p.nParticle()*p.mass()*magSqr(p.U() - u)
        );
    }
    uSqrAverage_->average(massAverage_);

    // sauter mean radius
    radiusAverage_() = volumeAverage_();
    weightAverage = 0;
    forAllConstIter(typename CloudType, cloud, iter)
    {
        const typename CloudType::parcelType& p = iter();
        const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), cloud.mesh());

        weightAverage.add
        (
            p.position(),
            tetIs,
            p.nParticle()*pow(p.volume(), 2.0/3.0)
        );
    }
    weightAverage.average();
    radiusAverage_->average(weightAverage);

    // collision frequency
    weightAverage = 0;
    forAllConstIter(typename CloudType, cloud, iter)
    {
        const typename CloudType::parcelType& p = iter();
        tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), cloud.mesh());

        const scalar a = volumeAverage_->interpolate(p.position(), tetIs);
        const scalar r = radiusAverage_->interpolate(p.position(), tetIs);
        const vector u = uAverage_->interpolate(p.position(), tetIs);

        const scalar f = 0.75*a/pow3(r)*sqr(0.5*p.d() + r)*mag(p.U() - u);

        frequencyAverage_->add(p.position(), tetIs, p.nParticle()*f*f);

        weightAverage.add(p.position(), tetIs, p.nParticle()*f);
    }
    frequencyAverage_->average(weightAverage);
}


template<class ParcelType>
template<class CloudType>
inline typename Foam::MPPICParcel<ParcelType>::template
TrackingData<CloudType>::trackPart
Foam::MPPICParcel<ParcelType>::TrackingData<CloudType>::part() const
{
    return part_;
}


template<class ParcelType>
template<class CloudType>
inline typename Foam::MPPICParcel<ParcelType>::template
TrackingData<CloudType>::trackPart&
Foam::MPPICParcel<ParcelType>::TrackingData<CloudType>::part()
{
    return part_;
}


// ************************************************************************* //
