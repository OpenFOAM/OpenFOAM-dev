/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "ConeInjection.H"
#include "Constant.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ConeInjection<CloudType>::setInjectionMethod()
{
    const word injectionMethod =
        this->coeffDict().template lookupOrDefault<word>
        (
            "injectionMethod",
            word::null
        );

    if (injectionMethod == "point" || injectionMethod == word::null)
    {
        injectionMethod_ = imPoint;

        topoChange();
    }
    else if (injectionMethod == "disc")
    {
        injectionMethod_ = imDisc;

        dInner_ =
            this->coeffDict().template lookup<scalar>("dInner", dimLength);
        dOuter_ =
            this->coeffDict().template lookup<scalar>("dOuter", dimLength);
    }
    else
    {
        FatalErrorInFunction
            << "injectionMethod must be either 'point' or 'disc'"
            << exit(FatalError);
    }
}


template<class CloudType>
void Foam::ConeInjection<CloudType>::setFlowType()
{
    const word flowType =
        this->coeffDict().template lookupOrDefault<word>
        (
            "flowType",
            word::null
        );

    if (flowType == "constantVelocity" || flowType == word::null)
    {
        flowType_ = ftConstantVelocity;

        Umag_.reset
        (
            Function1<scalar>::New
            (
                "Umag",
                this->owner().db().time().userUnits(),
                dimVelocity,
                this->coeffDict()
            ).ptr()
        );
    }
    else if (flowType == "pressureDrivenVelocity")
    {
        flowType_ = ftPressureDrivenVelocity;

        Pinj_.reset
        (
            Function1<scalar>::New
            (
                "Pinj",
                this->owner().db().time().userUnits(),
                dimPressure,
                this->coeffDict()
            ).ptr()
        );
    }
    else if (flowType == "flowRateAndDischarge")
    {
        flowType_ = ftFlowRateAndDischarge;

        dInner_ =
            this->coeffDict().template lookup<scalar>("dInner", dimLength);
        dOuter_ =
            this->coeffDict().template lookup<scalar>("dOuter", dimLength);

        Cd_.reset
        (
            Function1<scalar>::New
            (
                "Cd",
                this->owner().db().time().userUnits(),
                dimless,
                this->coeffDict()
            ).ptr()
        );
    }
    else
    {
        FatalErrorInFunction
            << "flowType must be either 'constantVelocity', "
            << "'pressureDrivenVelocity' or 'flowRateAndDischarge'"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConeInjection<CloudType>::ConeInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    injectionMethod_(imPoint),
    flowType_(ftConstantVelocity),
    position_
    (
        Function1<vector>::New
        (
            "position",
            this->owner().db().time().userUnits(),
            dimLength,
            this->coeffDict()
        )
    ),
    direction_
    (
        Function1<vector>::New
        (
            "direction",
            this->owner().db().time().userUnits(),
            dimless,
            this->coeffDict()
        )
    ),
    injectorCoordinates_(barycentric::uniform(NaN)),
    injectorCell_(-1),
    injectorTetFace_(-1),
    injectorTetPt_(-1),
    duration_(this->readDuration(dict, owner)),
    massFlowRate_(this->readMassFlowRate(dict, owner, duration_)),
    parcelsPerSecond_(this->readParcelsPerSecond(dict, owner)),
    thetaInner_
    (
        Function1<scalar>::New
        (
            "thetaInner",
            this->owner().db().time().userUnits(),
            unitDegrees,
            this->coeffDict()
        )
    ),
    thetaOuter_
    (
        Function1<scalar>::New
        (
            "thetaOuter",
            this->owner().db().time().userUnits(),
            unitDegrees,
            this->coeffDict()
        )
    ),
    sizeDistribution_
    (
        distribution::New
        (
            dimLength,
            this->coeffDict().subDict("sizeDistribution"),
            this->sizeSampleQ(),
            owner.rndGen().generator()
        )
    ),
    dInner_(vGreat),
    dOuter_(vGreat),
    Umag_(nullptr),
    Cd_(nullptr),
    Pinj_(nullptr)
{
    setInjectionMethod();

    setFlowType();

    topoChange();
}


template<class CloudType>
Foam::ConeInjection<CloudType>::ConeInjection
(
    const ConeInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    injectionMethod_(im.injectionMethod_),
    flowType_(im.flowType_),
    position_(im.position_, false),
    direction_(im.direction_, false),
    injectorCoordinates_(im.injectorCoordinates_),
    injectorCell_(im.injectorCell_),
    injectorTetFace_(im.injectorTetFace_),
    injectorTetPt_(im.injectorTetPt_),
    duration_(im.duration_),
    massFlowRate_(im.massFlowRate_, false),
    parcelsPerSecond_(im.parcelsPerSecond_, false),
    thetaInner_(im.thetaInner_, false),
    thetaOuter_(im.thetaOuter_, false),
    sizeDistribution_(im.sizeDistribution_, false),
    dInner_(im.dInner_),
    dOuter_(im.dOuter_),
    Umag_(im.Umag_, false),
    Cd_(im.Cd_, false),
    Pinj_(im.Pinj_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConeInjection<CloudType>::~ConeInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ConeInjection<CloudType>::topoChange()
{
    if (injectionMethod_ == imPoint && position_->constant())
    {
        vector position = position_->value(0);
        this->findCellAtPosition
        (
            position,
            injectorCoordinates_,
            injectorCell_,
            injectorTetFace_,
            injectorTetPt_
        );
    }
}


template<class CloudType>
Foam::scalar Foam::ConeInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::scalar Foam::ConeInjection<CloudType>::nParcelsToInject
(
    const scalar t0,
    const scalar t1
)
{
    if (t1 >= 0 && t0 < duration_)
    {
        return parcelsPerSecond_->integral(max(t0, 0), min(t1, duration_));
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::ConeInjection<CloudType>::massToInject
(
    const scalar t0,
    const scalar t1
)
{
    if (t1 >= 0 && t0 < duration_)
    {
        return massFlowRate_->integral(max(t0, 0), min(t1, duration_));
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
void Foam::ConeInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label,
    const scalar time,
    barycentric& coordinates,
    label& celli,
    label& tetFacei,
    label& tetPti,
    label& facei
)
{
    randomGenerator& rndGen = this->owner().rndGen();

    const scalar t = time - this->SOI_;

    switch (injectionMethod_)
    {
        case imPoint:
        {
            const point pos = position_->value(t);
            if (position_->constant())
            {
                coordinates = injectorCoordinates_;
                celli = injectorCell_;
                tetFacei = injectorTetFace_;
                tetPti = injectorTetPt_;
            }
            else
            {
                this->findCellAtPosition
                (
                    pos,
                    coordinates,
                    celli,
                    tetFacei,
                    tetPti,
                    false
                );
            }
            break;
        }
        case imDisc:
        {
            const scalar beta = twoPi*this->globalScalar01(rndGen);
            const scalar frac = this->globalScalar01(rndGen);
            const vector n = normalised(direction_->value(t));
            const vector t1 = normalised(perpendicular(n));
            const vector t2 = normalised(n ^ t1);
            const vector tanVec = t1*cos(beta) + t2*sin(beta);
            const scalar d = sqrt((1 - frac)*sqr(dInner_) + frac*sqr(dOuter_));
            const point pos = position_->value(t) + d/2*tanVec;
            this->findCellAtPosition
            (
                pos,
                coordinates,
                celli,
                tetFacei,
                tetPti,
                false
            );
            break;
        }
        default:
        {
            break;
        }
    }
}


template<class CloudType>
void Foam::ConeInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar time,
    typename CloudType::parcelType::trackingData& td,
    typename CloudType::parcelType& parcel
)
{
    const polyMesh& mesh = this->owner().mesh();

    randomGenerator& rndGen = this->owner().rndGen();

    const scalar t = time - this->SOI_;

    // Get the angle from the axis and the vector perpendicular from the axis.
    // If injecting at a point, then these are calculated from two new random
    // numbers. If a disc, then these calculations have already been done in
    // setPositionAndCell, so the angle and vector can be reverse engineered
    // from the position.
    scalar theta = vGreat;
    vector tanVec = vector::max;
    switch (injectionMethod_)
    {
        case imPoint:
        {
            const scalar beta = twoPi*rndGen.scalar01();
            const scalar frac = rndGen.scalar01();
            const vector n = normalised(direction_->value(t));
            const vector t1 = normalised(perpendicular(n));
            const vector t2 = normalised(n ^ t1);
            tanVec = t1*cos(beta) + t2*sin(beta);
            theta =
                sqrt
                (
                    (1 - frac)*sqr(thetaInner_->value(t))
                  + frac*sqr(thetaOuter_->value(t))
                );
            break;
        }
        case imDisc:
        {
            const scalar r = mag(parcel.position(mesh) - position_->value(t));
            const scalar frac = (2*r - dInner_)/(dOuter_ - dInner_);
            tanVec = normalised(parcel.position(mesh) - position_->value(t));
            theta =
                (1 - frac)*thetaInner_->value(t)
              + frac*thetaOuter_->value(t);
            break;
        }
        default:
        {
            break;
        }
    }

    // The direction of injection
    const vector dirVec =
        normalised
        (
            cos(theta)*normalised(direction_->value(t))
          + sin(theta)*tanVec
        );

    // Set the velocity
    switch (flowType_)
    {
        case ftConstantVelocity:
        {
            parcel.U() = Umag_->value(t)*dirVec;
            break;
        }
        case ftPressureDrivenVelocity:
        {
            const scalar pAmbient = this->owner().pAmbient();
            const scalar rho = parcel.rho();
            const scalar Umag = ::sqrt(2*(Pinj_->value(t) - pAmbient)/rho);
            parcel.U() = Umag*dirVec;
            break;
        }
        case ftFlowRateAndDischarge:
        {
            const scalar A = 0.25*pi*(sqr(dOuter_) - sqr(dInner_));
            const scalar Umag =
                massFlowRate_->value(t)/(parcel.rho()*Cd_->value(t)*A);
            parcel.U() = Umag*dirVec;
            break;
        }
        default:
        {
            break;
        }
    }

    // Set the particle diameter
    parcel.d() = sizeDistribution_->sample();
}


template<class CloudType>
bool Foam::ConeInjection<CloudType>::fullyDescribed() const
{
    return false;
}


// ************************************************************************* //
