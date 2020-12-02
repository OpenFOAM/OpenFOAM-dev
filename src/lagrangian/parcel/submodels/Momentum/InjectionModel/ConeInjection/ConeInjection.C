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

#include "ConeInjection.H"
#include "TimeFunction1.H"
#include "Constant.H"
#include "mathematicalConstants.H"
#include "unitConversion.H"

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

        updateMesh();
    }
    else if (injectionMethod == "disc")
    {
        injectionMethod_ = imDisc;

        this->coeffDict().lookup("dInner") >> dInner_;
        this->coeffDict().lookup("dOuter") >> dOuter_;
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

        Umag_.reset(this->coeffDict());
    }
    else if (flowType == "pressureDrivenVelocity")
    {
        flowType_ = ftPressureDrivenVelocity;

        Pinj_.reset(this->coeffDict());
    }
    else if (flowType == "flowRateAndDischarge")
    {
        flowType_ = ftFlowRateAndDischarge;

        this->coeffDict().lookup("dInner") >> dInner_;
        this->coeffDict().lookup("dOuter") >> dOuter_;

        Cd_.reset(this->coeffDict());
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
        TimeFunction1<vector>
        (
            owner.db().time(),
            "position",
            this->coeffDict()
        )
    ),
    positionIsConstant_(isA<Function1s::Constant<vector>>(position_)),
    direction_
    (
        TimeFunction1<vector>
        (
            owner.db().time(),
            "direction",
            this->coeffDict()
        )
    ),
    injectorCell_(-1),
    injectorTetFace_(-1),
    injectorTetPt_(-1),
    duration_(this->coeffDict().template lookup<scalar>("duration")),
    parcelsPerSecond_
    (
        this->coeffDict().template lookup<scalar>("parcelsPerSecond")
    ),
    flowRateProfile_
    (
        TimeFunction1<scalar>
        (
            owner.db().time(),
            "flowRateProfile",
            this->coeffDict()
        )
    ),
    thetaInner_
    (
        TimeFunction1<scalar>
        (
            owner.db().time(),
            "thetaInner",
            this->coeffDict()
        )
    ),
    thetaOuter_
    (
        TimeFunction1<scalar>
        (
            owner.db().time(),
            "thetaOuter",
            this->coeffDict()
        )
    ),
    sizeDistribution_
    (
        distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"), owner.rndGen()
        )
    ),
    dInner_(vGreat),
    dOuter_(vGreat),
    Umag_(owner.db().time(), "Umag"),
    Cd_(owner.db().time(), "Cd"),
    Pinj_(owner.db().time(), "Pinj")
{
    duration_ = owner.db().time().userTimeToTime(duration_);

    setInjectionMethod();

    setFlowType();

    // Set total volume to inject
    this->volumeTotal_ = flowRateProfile_.integral(0, duration_);

    updateMesh();
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
    position_(im.position_),
    positionIsConstant_(im.positionIsConstant_),
    direction_(im.direction_),
    injectorCell_(im.injectorCell_),
    injectorTetFace_(im.injectorTetFace_),
    injectorTetPt_(im.injectorTetPt_),
    duration_(im.duration_),
    parcelsPerSecond_(im.parcelsPerSecond_),
    flowRateProfile_(im.flowRateProfile_),
    thetaInner_(im.thetaInner_),
    thetaOuter_(im.thetaOuter_),
    sizeDistribution_(im.sizeDistribution_().clone().ptr()),
    dInner_(im.dInner_),
    dOuter_(im.dOuter_),
    Umag_(im.Umag_),
    Cd_(im.Cd_),
    Pinj_(im.Pinj_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConeInjection<CloudType>::~ConeInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ConeInjection<CloudType>::updateMesh()
{
    if (injectionMethod_ == imPoint && positionIsConstant_)
    {
        vector position = position_.value(0);
        this->findCellAtPosition
        (
            injectorCell_,
            injectorTetFace_,
            injectorTetPt_,
            position
        );
    }
}


template<class CloudType>
Foam::scalar Foam::ConeInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::ConeInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if (time0 >= 0 && time0 < duration_)
    {
        //// Standard calculation
        //return floor(parcelsPerSecond_*(time1 - time0));

        // Modified calculation to make numbers exact
        return floor(parcelsPerSecond_*time1 - this->parcelsAddedTotal());
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::ConeInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    if (time0 >= 0 && time0 < duration_)
    {
        return flowRateProfile_.integral(time0, time1);
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
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    Random& rndGen = this->owner().rndGen();

    const scalar t = time - this->SOI_;

    switch (injectionMethod_)
    {
        case imPoint:
        {
            position = position_.value(t);
            if (positionIsConstant_)
            {
                cellOwner = injectorCell_;
                tetFacei = injectorTetFace_;
                tetPti = injectorTetPt_;
            }
            else
            {
                this->findCellAtPosition
                (
                    cellOwner,
                    tetFacei,
                    tetPti,
                    position,
                    false
                );
            }
            break;
        }
        case imDisc:
        {
            const scalar beta = twoPi*rndGen.globalScalar01();
            const scalar frac = rndGen.globalScalar01();
            const vector n = normalised(direction_.value(t));
            const vector t1 = normalised(perpendicular(n));
            const vector t2 = normalised(n ^ t1);
            const vector tanVec = t1*cos(beta) + t2*sin(beta);
            const scalar d = sqrt((1 - frac)*sqr(dInner_) + frac*sqr(dOuter_));
            position = position_.value(t) + d/2*tanVec;
            this->findCellAtPosition
            (
                cellOwner,
                tetFacei,
                tetPti,
                position,
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
    typename CloudType::parcelType& parcel
)
{
    Random& rndGen = this->owner().rndGen();

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
            const vector n = normalised(direction_.value(t));
            const vector t1 = normalised(perpendicular(n));
            const vector t2 = normalised(n ^ t1);
            tanVec = t1*cos(beta) + t2*sin(beta);
            theta =
                degToRad
                (
                    sqrt
                    (
                        (1 - frac)*sqr(thetaInner_.value(t))
                        + frac*sqr(thetaOuter_.value(t))
                    )
                );
            break;
        }
        case imDisc:
        {
            const scalar r = mag(parcel.position() - position_.value(t));
            const scalar frac = (2*r - dInner_)/(dOuter_ - dInner_);
            tanVec = normalised(parcel.position() - position_.value(t));
            theta =
                degToRad
                (
                    (1 - frac)*thetaInner_.value(t)
                    + frac*thetaOuter_.value(t)
                );
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
            cos(theta)*normalised(direction_.value(t))
          + sin(theta)*tanVec
        );

    // Set the velocity
    switch (flowType_)
    {
        case ftConstantVelocity:
        {
            parcel.U() = Umag_.value(t)*dirVec;
            break;
        }
        case ftPressureDrivenVelocity:
        {
            const scalar pAmbient = this->owner().pAmbient();
            const scalar rho = parcel.rho();
            const scalar Umag = ::sqrt(2*(Pinj_.value(t) - pAmbient)/rho);
            parcel.U() = Umag*dirVec;
            break;
        }
        case ftFlowRateAndDischarge:
        {
            const scalar A = 0.25*pi*(sqr(dOuter_) - sqr(dInner_));
            const scalar massFlowRate =
                this->massTotal()*flowRateProfile_.value(t)/this->volumeTotal();
            const scalar Umag =
                massFlowRate/(parcel.rho()*Cd_.value(t)*A);
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


template<class CloudType>
bool Foam::ConeInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
