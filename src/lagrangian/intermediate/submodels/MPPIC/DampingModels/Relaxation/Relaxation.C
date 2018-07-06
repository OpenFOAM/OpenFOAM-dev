/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "Relaxation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DampingModels::Relaxation<CloudType>::Relaxation
(
    const dictionary& dict,
    CloudType& owner
)
:
    DampingModel<CloudType>(dict, owner, typeName),
    uAverage_(nullptr),
    oneByTimeScaleAverage_(nullptr)
{}


template<class CloudType>
Foam::DampingModels::Relaxation<CloudType>::Relaxation
(
    const Relaxation<CloudType>& cm
)
:
    DampingModel<CloudType>(cm),
    uAverage_(nullptr),
    oneByTimeScaleAverage_(cm.oneByTimeScaleAverage_->clone())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DampingModels::Relaxation<CloudType>::
~Relaxation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::DampingModels::Relaxation<CloudType>::cacheFields(const bool store)
{
    if (store)
    {
        const fvMesh& mesh = this->owner().mesh();
        const word& cloudName = this->owner().name();

        const AveragingMethod<scalar>& volumeAverage =
            mesh.lookupObject<AveragingMethod<scalar>>
            (
                cloudName + ":volumeAverage"
            );
        const AveragingMethod<scalar>& radiusAverage =
            mesh.lookupObject<AveragingMethod<scalar>>
            (
                cloudName + ":radiusAverage"
            );
        const AveragingMethod<vector>& uAverage =
            mesh.lookupObject<AveragingMethod<vector>>
            (
                cloudName + ":uAverage"
            );
        const AveragingMethod<scalar>& uSqrAverage =
            mesh.lookupObject<AveragingMethod<scalar>>
            (
                cloudName + ":uSqrAverage"
            );
        const AveragingMethod<scalar>& frequencyAverage =
            mesh.lookupObject<AveragingMethod<scalar>>
            (
                cloudName + ":frequencyAverage"
            );

        uAverage_ = &uAverage;

        oneByTimeScaleAverage_.reset
        (
            AveragingMethod<scalar>::New
            (
                IOobject
                (
                    cloudName + ":oneByTimeScaleAverage",
                    this->owner().db().time().timeName(),
                    mesh
                ),
                this->owner().solution().dict(),
                mesh
            ).ptr()
        );

        oneByTimeScaleAverage_() =
        (
            this->timeScaleModel_->oneByTau
            (
                volumeAverage,
                radiusAverage,
                uSqrAverage,
                frequencyAverage
            )
        )();
    }
    else
    {
        uAverage_ = nullptr;
        oneByTimeScaleAverage_.clear();
    }
}


template<class CloudType>
Foam::vector Foam::DampingModels::Relaxation<CloudType>::velocityCorrection
(
    typename CloudType::parcelType& p,
    const scalar deltaT
) const
{
    const tetIndices tetIs(p.currentTetIndices());

    const scalar x =
        deltaT*oneByTimeScaleAverage_->interpolate(p.coordinates(), tetIs);

    const vector u = uAverage_->interpolate(p.coordinates(), tetIs);

    return (u - p.U())*x/(x + 2.0);
}


// ************************************************************************* //
