/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "Explicit.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PackingModels::Explicit<CloudType>::Explicit
(
    const dictionary& dict,
    CloudType& owner
)
:
    PackingModel<CloudType>(dict, owner, typeName),
    stressAverage_(NULL),
    correctionLimiting_
    (
        CorrectionLimitingMethod::New
        (
            this->coeffDict().subDict(CorrectionLimitingMethod::typeName)
        )
    )
{}


template<class CloudType>
Foam::PackingModels::Explicit<CloudType>::Explicit
(
    const Explicit<CloudType>& cm
)
:
    PackingModel<CloudType>(cm),
    stressAverage_(cm.stressAverage_->clone()),
    correctionLimiting_
    (
        cm.correctionLimiting_->clone()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PackingModels::Explicit<CloudType>::~Explicit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PackingModels::Explicit<CloudType>::cacheFields(const bool store)
{
    PackingModel<CloudType>::cacheFields(store);

    if (store)
    {
        const fvMesh& mesh = this->owner().mesh();
        const word& cloudName = this->owner().name();

        const AveragingMethod<scalar>& volumeAverage =
            mesh.lookupObject<AveragingMethod<scalar> >
            (
                cloudName + ":volumeAverage"
            );
        const AveragingMethod<scalar>& rhoAverage =
            mesh.lookupObject<AveragingMethod<scalar> >
            (
                cloudName + ":rhoAverage"
            );
        const AveragingMethod<vector>& uAverage =
            mesh.lookupObject<AveragingMethod<vector> >
            (
                cloudName + ":uAverage"
            );
        const AveragingMethod<scalar>& uSqrAverage =
            mesh.lookupObject<AveragingMethod<scalar> >
            (
                cloudName + ":uSqrAverage"
            );

        volumeAverage_ = &volumeAverage;
        uAverage_ = &uAverage;

        stressAverage_.reset
        (
            AveragingMethod<scalar>::New
            (
                IOobject
                (
                    cloudName + ":stressAverage",
                    this->owner().db().time().timeName(),
                    mesh
                ),
                this->owner().solution().dict(),
                mesh
            ).ptr()
        );

        stressAverage_() =
            this->particleStressModel_->tau
            (
                *volumeAverage_,
                rhoAverage,
                uSqrAverage
            )();
    }
    else
    {
        volumeAverage_ = NULL;
        uAverage_ = NULL;
        stressAverage_.clear();
    }
}


template<class CloudType>
Foam::vector Foam::PackingModels::Explicit<CloudType>::velocityCorrection
(
    typename CloudType::parcelType& p,
    const scalar deltaT
) const
{
    const fvMesh& mesh = this->owner().mesh();
    const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), mesh);

    // interpolated quantities
    const scalar alpha =
        this->volumeAverage_->interpolate(p.position(), tetIs);
    const vector alphaGrad =
        this->volumeAverage_->interpolateGrad(p.position(), tetIs);
    const vector uMean =
        this->uAverage_->interpolate(p.position(), tetIs);

    // stress gradient
    const vector tauGrad =
        stressAverage_->interpolateGrad(p.position(), tetIs);

    // parcel relative velocity
    const vector uRelative = p.U() - uMean;

    // correction velocity
    vector dU = vector::zero;

    //// existing forces
    //const scalar Re = p.Re(p.U(), p.d(), p.rhoc(), p.muc());
    //const typename CloudType::forceType& forces = this->owner().forces();
    //const forceSuSp F =
    //    forces.calcCoupled(p, deltaT, p.mass(), Re, p.muc())
    //  + forces.calcNonCoupled(p, deltaT, p.mass(), Re, p.muc());

    // correction velocity
    if ((uRelative & alphaGrad) > 0)
    {
        dU = - deltaT*tauGrad/(p.rho()*alpha/* + deltaT*F.Sp()*/);
    }

    // apply the velocity limiters
    return
        correctionLimiting_->limitedVelocity
        (
            p.U(),
            dU,
            uMean
        );
}


// ************************************************************************* //
