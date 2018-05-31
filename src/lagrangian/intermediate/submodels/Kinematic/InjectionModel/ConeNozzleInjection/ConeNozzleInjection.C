/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "ConeNozzleInjection.H"
#include "TimeFunction1.H"
#include "mathematicalConstants.H"
#include "distributionModel.H"

using namespace Foam::constant;

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ConeNozzleInjection<CloudType>::setInjectionMethod()
{
    word injectionMethodType = this->coeffDict().lookup("injectionMethod");
    if (injectionMethodType == "disc")
    {
        injectionMethod_ = imDisc;
    }
    else if (injectionMethodType == "point")
    {
        injectionMethod_ = imPoint;

        // Set/cache the injector cell
        this->findCellAtPosition
        (
            injectorCell_,
            tetFacei_,
            tetPti_,
            position_,
            false
        );
    }
    else
    {
        FatalErrorInFunction
            << "injectionMethod must be either 'point' or 'disc'"
            << exit(FatalError);
    }
}


template<class CloudType>
void Foam::ConeNozzleInjection<CloudType>::setFlowType()
{
    word flowType = this->coeffDict().lookup("flowType");
    if (flowType == "constantVelocity")
    {
        this->coeffDict().lookup("UMag") >> UMag_;
        flowType_ = ftConstantVelocity;
    }
    else if (flowType == "pressureDrivenVelocity")
    {
        Pinj_.reset(this->coeffDict());
        flowType_ = ftPressureDrivenVelocity;
    }
    else if (flowType == "flowRateAndDischarge")
    {
        Cd_.reset(this->coeffDict());
        flowType_ = ftFlowRateAndDischarge;
    }
    else
    {
        FatalErrorInFunction
            << "flowType must be either 'constantVelocity', "
            <<"'pressureDrivenVelocity' or 'flowRateAndDischarge'"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConeNozzleInjection<CloudType>::ConeNozzleInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    injectionMethod_(imPoint),
    flowType_(ftConstantVelocity),
    outerDiameter_(readScalar(this->coeffDict().lookup("outerDiameter"))),
    innerDiameter_(readScalar(this->coeffDict().lookup("innerDiameter"))),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    position_(this->coeffDict().lookup("position")),
    injectorCell_(-1),
    tetFacei_(-1),
    tetPti_(-1),
    direction_(this->coeffDict().lookup("direction")),
    parcelsPerSecond_
    (
        readScalar(this->coeffDict().lookup("parcelsPerSecond"))
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
            this->coeffDict().subDict("sizeDistribution"),
            owner.rndGen()
        )
    ),
    tanVec1_(Zero),
    tanVec2_(Zero),
    normal_(Zero),

    UMag_(0.0),
    Cd_(owner.db().time(), "Cd"),
    Pinj_(owner.db().time(), "Pinj")
{
    if (innerDiameter_ >= outerDiameter_)
    {
        FatalErrorInFunction
         << "innerNozzleDiameter >= outerNozzleDiameter" << nl
         << exit(FatalError);
    }

    duration_ = owner.db().time().userTimeToTime(duration_);

    setInjectionMethod();

    setFlowType();

    // Normalise direction vector
    direction_ /= mag(direction_);

    // Determine direction vectors tangential to direction
    tanVec1_ = normalised(perpendicular(direction_));
    tanVec2_ = normalised(direction_ ^ tanVec1_);

    // Set total volume to inject
    this->volumeTotal_ = flowRateProfile_.integrate(0.0, duration_);

    updateMesh();
}


template<class CloudType>
Foam::ConeNozzleInjection<CloudType>::ConeNozzleInjection
(
    const ConeNozzleInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    injectionMethod_(im.injectionMethod_),
    flowType_(im.flowType_),
    outerDiameter_(im.outerDiameter_),
    innerDiameter_(im.innerDiameter_),
    duration_(im.duration_),
    position_(im.position_),
    injectorCell_(im.injectorCell_),
    tetFacei_(im.tetFacei_),
    tetPti_(im.tetPti_),
    direction_(im.direction_),
    parcelsPerSecond_(im.parcelsPerSecond_),
    flowRateProfile_(im.flowRateProfile_),
    thetaInner_(im.thetaInner_),
    thetaOuter_(im.thetaOuter_),
    sizeDistribution_(im.sizeDistribution_().clone().ptr()),
    tanVec1_(im.tanVec1_),
    tanVec2_(im.tanVec2_),
    normal_(im.normal_),
    UMag_(im.UMag_),
    Cd_(im.Cd_),
    Pinj_(im.Pinj_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConeNozzleInjection<CloudType>::~ConeNozzleInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ConeNozzleInjection<CloudType>::updateMesh()
{
    // Set/cache the injector cells
    switch (injectionMethod_)
    {
        case imPoint:
        {
            this->findCellAtPosition
            (
                injectorCell_,
                tetFacei_,
                tetPti_,
                position_
            );
        }
        default:
        {}
    }
}


template<class CloudType>
Foam::scalar Foam::ConeNozzleInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::ConeNozzleInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return floor((time1 - time0)*parcelsPerSecond_);
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::ConeNozzleInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return flowRateProfile_.integrate(time0, time1);
    }
    else
    {
        return 0.0;
    }
}


template<class CloudType>
void Foam::ConeNozzleInjection<CloudType>::setPositionAndCell
(
    const label,
    const label,
    const scalar,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    cachedRandom& rndGen = this->owner().rndGen();

    scalar beta = mathematical::twoPi*rndGen.globalSample01<scalar>();
    normal_ = tanVec1_*cos(beta) + tanVec2_*sin(beta);

    switch (injectionMethod_)
    {
        case imPoint:
        {
            position = position_;
            cellOwner = injectorCell_;
            tetFacei = tetFacei_;
            tetPti = tetPti_;

            break;
        }
        case imDisc:
        {
            scalar frac = rndGen.globalSample01<scalar>();
            scalar dr = outerDiameter_ - innerDiameter_;
            scalar r = 0.5*(innerDiameter_ + frac*dr);
            position = position_ + r*normal_;

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
            FatalErrorInFunction
             << "Unknown injectionMethod type" << nl
             << exit(FatalError);
        }
    }
}


template<class CloudType>
void Foam::ConeNozzleInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar time,
    typename CloudType::parcelType& parcel
)
{
    cachedRandom& rndGen = this->owner().rndGen();

    // set particle velocity
    const scalar deg2Rad = mathematical::pi/180.0;

    scalar t = time - this->SOI_;
    scalar ti = thetaInner_.value(t);
    scalar to = thetaOuter_.value(t);
    scalar coneAngle = rndGen.sample01<scalar>()*(to - ti) + ti;

    coneAngle *= deg2Rad;
    scalar alpha = sin(coneAngle);
    scalar dcorr = cos(coneAngle);

    vector normal = alpha*normal_;
    vector dirVec = dcorr*direction_;
    dirVec += normal;
    dirVec /= mag(dirVec);

    switch (flowType_)
    {
        case ftConstantVelocity:
        {
            parcel.U() = UMag_*dirVec;
            break;
        }
        case ftPressureDrivenVelocity:
        {
            scalar pAmbient = this->owner().pAmbient();
            scalar rho = parcel.rho();
            scalar UMag = ::sqrt(2.0*(Pinj_.value(t) - pAmbient)/rho);
            parcel.U() = UMag*dirVec;
            break;
        }
        case ftFlowRateAndDischarge:
        {
            scalar Ao = 0.25*mathematical::pi*outerDiameter_*outerDiameter_;
            scalar Ai = 0.25*mathematical::pi*innerDiameter_*innerDiameter_;
            scalar massFlowRate =
                this->massTotal()
               *flowRateProfile_.value(t)
               /this->volumeTotal();

            scalar Umag = massFlowRate/(parcel.rho()*Cd_.value(t)*(Ao - Ai));
            parcel.U() = Umag*dirVec;
            break;
        }
        default:
        {
        }
    }

    // set particle diameter
    parcel.d() = sizeDistribution_->sample();
}


template<class CloudType>
bool Foam::ConeNozzleInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::ConeNozzleInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
