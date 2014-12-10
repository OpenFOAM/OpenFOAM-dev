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

#include "ThermoSurfaceFilm.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "Pstream.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
Foam::wordList Foam::ThermoSurfaceFilm<CloudType>::interactionTypeNames_
(
    IStringStream
    (
        "(absorb bounce splashBai)"
    )()
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class CloudType>
typename Foam::ThermoSurfaceFilm<CloudType>::interactionType
Foam::ThermoSurfaceFilm<CloudType>::interactionTypeEnum(const word& it) const
{
    forAll(interactionTypeNames_, i)
    {
        if (interactionTypeNames_[i] == it)
        {
            return interactionType(i);
        }
    }

    FatalErrorIn
    (
        "ThermoSurfaceFilm<CloudType>::interactionType "
        "ThermoSurfaceFilm<CloudType>::interactionTypeEnum"
        "("
            "const word& it"
        ") const"
    )   << "Unknown interaction type " << it
        << ". Valid interaction types include: " << interactionTypeNames_
        << abort(FatalError);

    return interactionType(0);
}


template<class CloudType>
Foam::word Foam::ThermoSurfaceFilm<CloudType>::interactionTypeStr
(
    const interactionType& it
) const
{
    if (it >= interactionTypeNames_.size())
    {
        FatalErrorIn
        (
            "ThermoSurfaceFilm<CloudType>::interactionType "
            "ThermoSurfaceFilm<CloudType>::interactionTypeStr"
            "("
                "const interactionType& it"
            ") const"
        )   << "Unknown interaction type enumeration" << abort(FatalError);
    }

    return interactionTypeNames_[it];
}


template<class CloudType>
Foam::vector Foam::ThermoSurfaceFilm<CloudType>::tangentVector
(
    const vector& v
) const
{
    vector tangent = vector::zero;
    scalar magTangent = 0.0;

    while (magTangent < SMALL)
    {
        vector vTest = rndGen_.sample01<vector>();
        tangent = vTest - (vTest & v)*v;
        magTangent = mag(tangent);
    }

    return tangent/magTangent;
}


template<class CloudType>
Foam::vector Foam::ThermoSurfaceFilm<CloudType>::splashDirection
(
    const vector& tanVec1,
    const vector& tanVec2,
    const vector& nf
) const
{
    // azimuthal angle [rad]
    const scalar phiSi = twoPi*rndGen_.sample01<scalar>();

    // ejection angle [rad]
    const scalar thetaSi = pi/180.0*(rndGen_.sample01<scalar>()*(50 - 5) + 5);

    // direction vector of new parcel
    const scalar alpha = sin(thetaSi);
    const scalar dcorr = cos(thetaSi);
    const vector normal = alpha*(tanVec1*cos(phiSi) + tanVec2*sin(phiSi));
    vector dirVec = dcorr*nf;
    dirVec += normal;

    return dirVec/mag(dirVec);
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::absorbInteraction
(
    regionModels::surfaceFilmModels::surfaceFilmModel& filmModel,
    const parcelType& p,
    const polyPatch& pp,
    const label faceI,
    const scalar mass,
    bool& keepParticle
)
{
    if (debug)
    {
        Info<< "Parcel " << p.origId() << " absorbInteraction" << endl;
    }

    // Patch face normal
    const vector& nf = pp.faceNormals()[faceI];

    // Patch velocity
    const vector& Up = this->owner().U().boundaryField()[pp.index()][faceI];

    // Relative parcel velocity
    const vector Urel = p.U() - Up;

    // Parcel normal velocity
    const vector Un = nf*(Urel & nf);

    // Parcel tangential velocity
    const vector Ut = Urel - Un;

    filmModel.addSources
    (
        pp.index(),
        faceI,
        mass,                           // mass
        mass*Ut,                        // tangential momentum
        mass*mag(Un),                   // impingement pressure
        mass*p.hs()                     // energy
    );

    this->nParcelsTransferred()++;

    keepParticle = false;
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::bounceInteraction
(
    parcelType& p,
    const polyPatch& pp,
    const label faceI,
    bool& keepParticle
) const
{
    if (debug)
    {
        Info<< "Parcel " << p.origId() << " bounceInteraction" << endl;
    }

    // Patch face normal
    const vector& nf = pp.faceNormals()[faceI];

    // Patch velocity
    const vector& Up = this->owner().U().boundaryField()[pp.index()][faceI];

    // Relative parcel velocity
    const vector Urel = p.U() - Up;

    // Flip parcel normal velocity component
    p.U() -= 2.0*nf*(Urel & nf);

    keepParticle = true;
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::drySplashInteraction
(
    regionModels::surfaceFilmModels::surfaceFilmModel& filmModel,
    const parcelType& p,
    const polyPatch& pp,
    const label faceI,
    bool& keepParticle
)
{
    if (debug)
    {
        Info<< "Parcel " << p.origId() << " drySplashInteraction" << endl;
    }

    const liquidProperties& liq = thermo_.liquids().properties()[0];

    // Patch face velocity and normal
    const vector& Up = this->owner().U().boundaryField()[pp.index()][faceI];
    const vector& nf = pp.faceNormals()[faceI];

    // local pressure
    const scalar pc = thermo_.thermo().p()[p.cell()];

    // Retrieve parcel properties
    const scalar m = p.mass()*p.nParticle();
    const scalar rho = p.rho();
    const scalar d = p.d();
    const scalar sigma = liq.sigma(pc, p.T());
    const scalar mu = liq.mu(pc, p.T());
    const vector Urel = p.U() - Up;
    const vector Un = nf*(Urel & nf);

    // Laplace number
    const scalar La = rho*sigma*d/sqr(mu);

    // Weber number
    const scalar We = rho*magSqr(Un)*d/sigma;

    // Critical Weber number
    const scalar Wec = Adry_*pow(La, -0.183);

    if (We < Wec) // adhesion - assume absorb
    {
        absorbInteraction(filmModel, p, pp, faceI, m, keepParticle);
    }
    else // splash
    {
        // ratio of incident mass to splashing mass
        const scalar mRatio = 0.2 + 0.6*rndGen_.sample01<scalar>();
        splashInteraction
            (filmModel, p, pp, faceI, mRatio, We, Wec, sigma, keepParticle);
    }
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::wetSplashInteraction
(
    regionModels::surfaceFilmModels::surfaceFilmModel& filmModel,
    parcelType& p,
    const polyPatch& pp,
    const label faceI,
    bool& keepParticle
)
{
    if (debug)
    {
        Info<< "Parcel " << p.origId() << " wetSplashInteraction" << endl;
    }

    const liquidProperties& liq = thermo_.liquids().properties()[0];

    // Patch face velocity and normal
    const vector& Up = this->owner().U().boundaryField()[pp.index()][faceI];
    const vector& nf = pp.faceNormals()[faceI];

    // local pressure
    const scalar pc = thermo_.thermo().p()[p.cell()];

    // Retrieve parcel properties
    const scalar m = p.mass()*p.nParticle();
    const scalar rho = p.rho();
    const scalar d = p.d();
    vector& U = p.U();
    const scalar sigma = liq.sigma(pc, p.T());
    const scalar mu = liq.mu(pc, p.T());
    const vector Urel = p.U() - Up;
    const vector Un = nf*(Urel & nf);
    const vector Ut = Urel - Un;

    // Laplace number
    const scalar La = rho*sigma*d/sqr(mu);

    // Weber number
    const scalar We = rho*magSqr(Un)*d/sigma;

    // Critical Weber number
    const scalar Wec = Awet_*pow(La, -0.183);

    if (We < 1) // adhesion - assume absorb
    {
        absorbInteraction(filmModel, p, pp, faceI, m, keepParticle);
    }
    else if ((We >= 1) && (We < 20)) // bounce
    {
        // incident angle of impingement
        const scalar theta = pi/2 - acos(U/mag(U) & nf);

        // restitution coefficient
        const scalar epsilon = 0.993 - theta*(1.76 - theta*(1.56 - theta*0.49));

        // update parcel velocity
        U = -epsilon*(Un) + 5/7*(Ut);

        keepParticle = true;
        return;
    }
    else if ((We >= 20) && (We < Wec)) // spread - assume absorb
    {
        absorbInteraction(filmModel, p, pp, faceI, m, keepParticle);
    }
    else    // splash
    {
        // ratio of incident mass to splashing mass
        // splash mass can be > incident mass due to film entrainment
        const scalar mRatio = 0.2 + 0.9*rndGen_.sample01<scalar>();
        splashInteraction
            (filmModel, p, pp, faceI, mRatio, We, Wec, sigma, keepParticle);
    }
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::splashInteraction
(
    regionModels::surfaceFilmModels::surfaceFilmModel& filmModel,
    const parcelType& p,
    const polyPatch& pp,
    const label faceI,
    const scalar mRatio,
    const scalar We,
    const scalar Wec,
    const scalar sigma,
    bool& keepParticle
)
{
    // Patch face velocity and normal
    const fvMesh& mesh = this->owner().mesh();
    const vector& Up = this->owner().U().boundaryField()[pp.index()][faceI];
    const vector& nf = pp.faceNormals()[faceI];

    // Determine direction vectors tangential to patch normal
    const vector tanVec1 = tangentVector(nf);
    const vector tanVec2 = nf^tanVec1;

    // Retrieve parcel properties
    const scalar np = p.nParticle();
    const scalar m = p.mass()*np;
    const scalar d = p.d();
    const vector Urel = p.U() - Up;
    const vector Un = nf*(Urel & nf);
    const vector Ut = Urel - Un;
    const vector& posC = mesh.C()[p.cell()];
    const vector& posCf = mesh.Cf().boundaryField()[pp.index()][faceI];

    // total mass of (all) splashed parcels
    const scalar mSplash = m*mRatio;

    // number of splashed particles per incoming particle
    const scalar Ns = 5.0*(We/Wec - 1.0);

    // average diameter of splashed particles
    const scalar dBarSplash = 1/cbrt(6.0)*cbrt(mRatio/Ns)*d + ROOTVSMALL;

    // cumulative diameter splash distribution
    const scalar dMax = 0.9*cbrt(mRatio)*d;
    const scalar dMin = 0.1*dMax;
    const scalar K = exp(-dMin/dBarSplash) - exp(-dMax/dBarSplash);

    // surface energy of secondary parcels [J]
    scalar ESigmaSec = 0;

    // sample splash distribution to determine secondary parcel diameters
    scalarList dNew(parcelsPerSplash_);
    scalarList npNew(parcelsPerSplash_);
    forAll(dNew, i)
    {
        const scalar y = rndGen_.sample01<scalar>();
        dNew[i] = -dBarSplash*log(exp(-dMin/dBarSplash) - y*K);
        npNew[i] = mRatio*np*pow3(d)/pow3(dNew[i])/parcelsPerSplash_;
        ESigmaSec += npNew[i]*sigma*p.areaS(dNew[i]);
    }

    // incident kinetic energy [J]
    const scalar EKIn = 0.5*m*magSqr(Urel);

    // incident surface energy [J]
    const scalar ESigmaIn = np*sigma*p.areaS(d);

    // dissipative energy
    const scalar Ed = max(0.8*EKIn, np*Wec/12*pi*sigma*sqr(d));

    // total energy [J]
    const scalar EKs = EKIn + ESigmaIn - ESigmaSec - Ed;

    // switch to absorb if insufficient energy for splash
    if (EKs <= 0)
    {
        absorbInteraction(filmModel, p, pp, faceI, m, keepParticle);
        return;
    }

    // helper variables to calculate magUns0
    const scalar logD = log(d);
    const scalar coeff2 = log(dNew[0]) - logD + ROOTVSMALL;
    scalar coeff1 = 0.0;
    forAll(dNew, i)
    {
        coeff1 += sqr(log(dNew[i]) - logD);
    }

    // magnitude of the normal velocity of the first splashed parcel
    const scalar magUns0 =
        sqrt(2.0*parcelsPerSplash_*EKs/mSplash/(1.0 + coeff1/sqr(coeff2)));

    // Set splashed parcel properties
    forAll(dNew, i)
    {
        const vector dirVec = splashDirection(tanVec1, tanVec2, -nf);

        // Create a new parcel by copying source parcel
        parcelType* pPtr = new parcelType(p);

        pPtr->origId() = pPtr->getNewParticleID();

        pPtr->origProc() = Pstream::myProcNo();

        if (splashParcelType_ >= 0)
        {
            pPtr->typeId() = splashParcelType_;
        }

        // perturb new parcels towards the owner cell centre
        pPtr->position() += 0.5*rndGen_.sample01<scalar>()*(posC - posCf);

        pPtr->nParticle() = npNew[i];

        pPtr->d() = dNew[i];

        pPtr->U() = dirVec*(mag(Cf_*Ut) + magUns0*(log(dNew[i]) - logD)/coeff2);

        // Apply correction to velocity for 2-D cases
        meshTools::constrainDirection(mesh, mesh.solutionD(), pPtr->U());

        // Add the new parcel
        this->owner().addParticle(pPtr);

        nParcelsSplashed_++;
    }

    // transfer remaining part of parcel to film 0 - splashMass can be -ve
    // if entraining from the film
    const scalar mDash = m - mSplash;
    absorbInteraction(filmModel, p, pp, faceI, mDash, keepParticle);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoSurfaceFilm<CloudType>::ThermoSurfaceFilm
(
    const dictionary& dict,
    CloudType& owner
)
:
    SurfaceFilmModel<CloudType>(dict, owner, typeName),
    rndGen_(owner.rndGen()),
    thermo_
    (
        owner.db().objectRegistry::template lookupObject<SLGThermo>("SLGThermo")
    ),
    TFilmPatch_(0),
    CpFilmPatch_(0),
    interactionType_
    (
        interactionTypeEnum(this->coeffDict().lookup("interactionType"))
    ),
    deltaWet_(0.0),
    splashParcelType_(0),
    parcelsPerSplash_(0),
    Adry_(0.0),
    Awet_(0.0),
    Cf_(0.0),
    nParcelsSplashed_(0)
{
    Info<< "    Applying " << interactionTypeStr(interactionType_)
        << " interaction model" << endl;

    if (interactionType_ == itSplashBai)
    {
        this->coeffDict().lookup("deltaWet") >> deltaWet_;
        splashParcelType_ =
            this->coeffDict().lookupOrDefault("splashParcelType", -1);
        parcelsPerSplash_ =
            this->coeffDict().lookupOrDefault("parcelsPerSplash", 2);
        this->coeffDict().lookup("Adry") >> Adry_;
        this->coeffDict().lookup("Awet") >> Awet_;
        this->coeffDict().lookup("Cf") >> Cf_;
    }
}


template<class CloudType>
Foam::ThermoSurfaceFilm<CloudType>::ThermoSurfaceFilm
(
    const ThermoSurfaceFilm<CloudType>& sfm
)
:
    SurfaceFilmModel<CloudType>(sfm),
    rndGen_(sfm.rndGen_),
    thermo_(sfm.thermo_),
    TFilmPatch_(sfm.TFilmPatch_),
    CpFilmPatch_(sfm.CpFilmPatch_),
    interactionType_(sfm.interactionType_),
    deltaWet_(sfm.deltaWet_),
    splashParcelType_(sfm.splashParcelType_),
    parcelsPerSplash_(sfm.parcelsPerSplash_),
    Adry_(sfm.Adry_),
    Awet_(sfm.Awet_),
    Cf_(sfm.Cf_),
    nParcelsSplashed_(sfm.nParcelsSplashed_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoSurfaceFilm<CloudType>::~ThermoSurfaceFilm()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::ThermoSurfaceFilm<CloudType>::transferParcel
(
    parcelType& p,
    const polyPatch& pp,
    bool& keepParticle
)
{
    // Retrieve the film model from the owner database
    regionModels::surfaceFilmModels::surfaceFilmModel& filmModel =
        const_cast<regionModels::surfaceFilmModels::surfaceFilmModel&>
        (
            this->owner().db().time().objectRegistry::template
                lookupObject<regionModels::surfaceFilmModels::surfaceFilmModel>
                (
                    "surfaceFilmProperties"
                )
        );

    const label patchI = pp.index();

    if (filmModel.isRegionPatch(patchI))
    {
        const label faceI = pp.whichFace(p.face());

        switch (interactionType_)
        {
            case itBounce:
            {
                bounceInteraction(p, pp, faceI, keepParticle);

                break;
            }
            case itAbsorb:
            {
                const scalar m = p.nParticle()*p.mass();
                absorbInteraction(filmModel, p, pp, faceI, m, keepParticle);

                break;
            }
            case itSplashBai:
            {
                bool dry = this->deltaFilmPatch_[patchI][faceI] < deltaWet_;

                if (dry)
                {
                    drySplashInteraction(filmModel, p, pp, faceI, keepParticle);
                }
                else
                {
                    wetSplashInteraction(filmModel, p, pp, faceI, keepParticle);
                }

                break;
            }
            default:
            {
                FatalErrorIn
                (
                    "bool ThermoSurfaceFilm<CloudType>::transferParcel"
                    "("
                        "parcelType&, "
                        "const polyPatch&, "
                        "bool&"
                    ")"
                )   << "Unknown interaction type enumeration"
                    << abort(FatalError);
            }
        }

        // transfer parcel/parcel interactions complete
        return true;
    }

    // parcel not interacting with film
    return false;
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::cacheFilmFields
(
    const label filmPatchI,
    const label primaryPatchI,
    const regionModels::surfaceFilmModels::surfaceFilmModel& filmModel
)
{
    SurfaceFilmModel<CloudType>::cacheFilmFields
    (
        filmPatchI,
        primaryPatchI,
        filmModel
    );

    TFilmPatch_ = filmModel.Ts().boundaryField()[filmPatchI];
    filmModel.toPrimary(filmPatchI, TFilmPatch_);

    CpFilmPatch_ = filmModel.Cp().boundaryField()[filmPatchI];
    filmModel.toPrimary(filmPatchI, CpFilmPatch_);
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::setParcelProperties
(
    parcelType& p,
    const label filmFaceI
) const
{
    SurfaceFilmModel<CloudType>::setParcelProperties(p, filmFaceI);

    // Set parcel properties
    p.T() = TFilmPatch_[filmFaceI];
    p.Cp() = CpFilmPatch_[filmFaceI];
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::info(Ostream& os)
{
    SurfaceFilmModel<CloudType>::info(os);

    label nSplash0 = this->template getModelProperty<label>("nParcelsSplashed");
    label nSplashTotal =
        nSplash0 + returnReduce(nParcelsSplashed_, sumOp<label>());

    os  << "    New film splash parcels         = " << nSplashTotal << endl;

    if (this->outputTime())
    {
        this->setModelProperty("nParcelsSplashed", nSplashTotal);
        nParcelsSplashed_ = 0;
    }
}


// ************************************************************************* //
