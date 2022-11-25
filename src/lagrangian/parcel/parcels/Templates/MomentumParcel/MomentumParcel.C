/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "MomentumParcel.H"
#include "forceSuSp.H"
#include "integrationScheme.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::label Foam::MomentumParcel<ParcelType>::maxTrackAttempts = 1;


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::MomentumParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    tetIndices tetIs = this->currentTetIndices(td.mesh);

    td.rhoc() = td.rhoInterp().interpolate(this->coordinates(), tetIs);

    if (td.rhoc() < cloud.constProps().rhoMin())
    {
        if (debug)
        {
            WarningInFunction
                << "Limiting observed density in cell " << this->cell()
                << " to " << cloud.constProps().rhoMin() <<  nl << endl;
        }

        td.rhoc() = cloud.constProps().rhoMin();
    }

    td.Uc() = td.UInterp().interpolate(this->coordinates(), tetIs);

    td.muc() = td.muInterp().interpolate(this->coordinates(), tetIs);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::MomentumParcel<ParcelType>::calcDispersion
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    td.Uc() = cloud.dispersion().update
    (
        dt,
        this->cell(),
        U_,
        td.Uc(),
        UTurb_,
        tTurb_
    );
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::MomentumParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    td.Uc() += cloud.UTransRef()[this->cell()]/massCell(td);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::MomentumParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = nParticle_;
    const scalar mass0 = mass();

    // Reynolds number
    const scalar Re = this->Re(td);


    // Sources
    //~~~~~~~~

    // Explicit momentum source for particle
    vector Su = Zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = Zero;


    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    this->U_ =
        calcVelocity(cloud, td, dt, Re, td.muc(), mass0, Su, dUTrans, Spu);


    // Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (cloud.solution().coupled())
    {
        // Update momentum transfer
        cloud.UTransRef()[this->cell()] += np0*dUTrans;

        // Update momentum transfer coefficient
        cloud.UCoeffRef()[this->cell()] += np0*Spu;
    }
}


template<class ParcelType>
template<class TrackCloudType>
const Foam::vector Foam::MomentumParcel<ParcelType>::calcVelocity
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar Re,
    const scalar mu,
    const scalar mass,
    const vector& Su,
    vector& dUTrans,
    scalar& Spu
) const
{
    const typename TrackCloudType::parcelType& p =
        static_cast<const typename TrackCloudType::parcelType&>(*this);
    typename TrackCloudType::parcelType::trackingData& ttd =
        static_cast<typename TrackCloudType::parcelType::trackingData&>(td);

    const typename TrackCloudType::forceType& forces = cloud.forces();

    // Momentum source due to particle forces
    const forceSuSp Fcp = forces.calcCoupled(p, ttd, dt, mass, Re, mu);
    const forceSuSp Fncp = forces.calcNonCoupled(p, ttd, dt, mass, Re, mu);
    const scalar massEff = forces.massEff(p, ttd, mass);

    /*
    // Proper splitting ...
    // Calculate the integration coefficients
    const vector acp = (Fcp.Sp()*td.Uc() + Fcp.Su())/massEff;
    const vector ancp = (Fncp.Sp()*td.Uc() + Fncp.Su() + Su)/massEff;
    const scalar bcp = Fcp.Sp()/massEff;
    const scalar bncp = Fncp.Sp()/massEff;

    // Integrate to find the new parcel velocity
    const vector deltaUcp =
        cloud.UIntegrator().partialDelta
        (
            U_, dt, acp + ancp, bcp + bncp, acp, bcp
        );
    const vector deltaUncp =
        cloud.UIntegrator().partialDelta
        (
            U_, dt, acp + ancp, bcp + bncp, ancp, bncp
        );
    const vector deltaT = deltaUcp + deltaUncp;
    */

    // Shortcut splitting assuming no implicit non-coupled force ...
    // Calculate the integration coefficients
    const vector acp = (Fcp.Sp()*td.Uc() + Fcp.Su())/massEff;
    const vector ancp = (Fncp.Su() + Su)/massEff;
    const scalar bcp = Fcp.Sp()/massEff;

    // Integrate to find the new parcel velocity
    const vector deltaU = cloud.UIntegrator().delta(U_, dt, acp + ancp, bcp);
    const vector deltaUncp = ancp*dt;
    const vector deltaUcp = deltaU - deltaUncp;

    // Calculate the new velocity and the momentum transfer terms
    vector Unew = U_ + deltaU;

    dUTrans -= massEff*deltaUcp;

    Spu = dt*Fcp.Sp();

    // Apply correction to velocity and dUTrans for reduced-D cases
    const polyMesh& mesh = cloud.pMesh();
    meshTools::constrainDirection(mesh, mesh.solutionD(), Unew);
    meshTools::constrainDirection(mesh, mesh.solutionD(), dUTrans);

    return Unew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::MomentumParcel<ParcelType>::MomentumParcel
(
    const MomentumParcel<ParcelType>& p
)
:
    ParcelType(p),
    moving_(p.moving_),
    typeId_(p.typeId_),
    nParticle_(p.nParticle_),
    d_(p.d_),
    dTarget_(p.dTarget_),
    U_(p.U_),
    rho_(p.rho_),
    age_(p.age_),
    tTurb_(p.tTurb_),
    UTurb_(p.UTurb_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
bool Foam::MomentumParcel<ParcelType>::move
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);
    typename TrackCloudType::parcelType::trackingData& ttd =
        static_cast<typename TrackCloudType::parcelType::trackingData&>(td);

    ttd.keepParticle = true;
    td.sendToProc = -1;

    const scalarField& cellLengthScale = cloud.cellLengthScale();
    const scalar maxCo = cloud.solution().maxCo();

    while
    (
        ttd.keepParticle
     && ttd.sendToProc == -1
     && p.stepFraction() < ttd.stepFractionRange().second()
    )
    {
        if (p.moving() && p.onFace())
        {
            cloud.functions().postFace(p);
        }

        // Cache the current position, cell and step-fraction
        const point start = p.position(td.mesh);
        const scalar sfrac = p.stepFraction();

        // Total displacement over the time-step
        const vector s = ttd.trackTime()*U_;

        // Cell length scale
        const scalar l = cellLengthScale[p.cell()];

        // Deviation from the mesh centre for reduced-D cases
        const vector d = p.deviationFromMeshCentre(td.mesh);

        // Fraction of the displacement to track in this loop. This is limited
        // to ensure that the both the time and distance tracked is less than
        // maxCo times the total value.
        scalar f = ttd.stepFractionRange().second() - p.stepFraction();
        f = min(f, maxCo);
        f = min(f, maxCo/min(max(mag(s)/l, rootSmall), rootGreat));
        if (p.moving())
        {
            // Track to the next face
            p.trackToFace(td.mesh, f*s - d, f);
        }
        else
        {
            // At present the only thing that sets moving_ to false is a stick
            // wall interaction. We want the position of the particle to remain
            // the same relative to the face that it is on. The local
            // coordinates therefore do not change. We still advance in time and
            // perform the relevant interactions with the fixed particle.
            p.stepFraction() += f;
        }

        const scalar dt = (p.stepFraction() - sfrac)*ttd.trackTime();

        // Avoid problems with extremely small timesteps
        if (dt > rootVSmall)
        {
            // Update cell based properties
            p.setCellValues(cloud, ttd);

            p.calcDispersion(cloud, ttd, dt);

            if (cloud.solution().cellValueSourceCorrection())
            {
                p.cellValueSourceCorrection(cloud, ttd, dt);
            }

            p.calc(cloud, ttd, dt);
        }

        p.age() += dt;

        cloud.functions().postMove(p, dt, start, ttd.keepParticle);

        if (p.moving() && p.onFace() && ttd.keepParticle)
        {
            cloud.functions().preFace(p);

            p.hitFace(f*s - d, f, cloud, ttd);
        }
    }

    return ttd.keepParticle;
}


template<class ParcelType>
void Foam::MomentumParcel<ParcelType>::transformProperties
(
    const transformer& transform
)
{
    ParcelType::transformProperties(transform);

    U_ = transform.transform(U_);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::MomentumParcel<ParcelType>::correctAfterParallelTransfer
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::correctAfterParallelTransfer(cloud, td);

    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);

    const polyPatch& pp = td.mesh.boundaryMesh()[td.sendToPatch];

    cloud.functions().postPatch(p, pp);
}


template<class ParcelType>
template<class TrackCloudType>
bool Foam::MomentumParcel<ParcelType>::hitPatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);

    const polyPatch& pp = td.mesh.boundaryMesh()[p.patch(td.mesh)];

    // Allow a surface film model to consume the parcel
    if (cloud.surfaceFilm().transferParcel(p, pp, td.keepParticle))
    {
        cloud.functions().postPatch(p, pp);
        return true;
    }

    // Pass to the patch interaction model
    if (cloud.patchInteraction().correct(p, pp, td.keepParticle))
    {
        cloud.functions().postPatch(p, pp);
        return true;
    }

    return false;
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::MomentumParcel<ParcelType>::hitWedgePatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::hitWedgePatch(cloud, td);

    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);

    const polyPatch& pp = td.mesh.boundaryMesh()[p.patch(td.mesh)];

    cloud.functions().postPatch(p, pp);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::MomentumParcel<ParcelType>::hitSymmetryPlanePatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::hitSymmetryPlanePatch(cloud, td);

    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);

    const polyPatch& pp = td.mesh.boundaryMesh()[p.patch(td.mesh)];

    cloud.functions().postPatch(p, pp);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::MomentumParcel<ParcelType>::hitSymmetryPatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::hitSymmetryPatch(cloud, td);

    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);

    const polyPatch& pp = td.mesh.boundaryMesh()[p.patch(td.mesh)];

    cloud.functions().postPatch(p, pp);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::MomentumParcel<ParcelType>::hitCyclicPatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::hitCyclicPatch(cloud, td);

    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);

    const polyPatch& pp = td.mesh.boundaryMesh()[p.patch(td.mesh)];

    cloud.functions().postPatch(p, pp);
}


template<class ParcelType>
template<class TrackCloudType>
bool Foam::MomentumParcel<ParcelType>::hitNonConformalCyclicPatch
(
    const vector& displacement,
    const scalar fraction,
    const label patchi,
    TrackCloudType& cloud,
    trackingData& td
)
{
    const bool result =
        ParcelType::hitNonConformalCyclicPatch
        (
            displacement,
            fraction,
            patchi,
            cloud,
            td
        );

    if (td.sendToProc == -1 && result)
    {
        typename TrackCloudType::parcelType& p =
            static_cast<typename TrackCloudType::parcelType&>(*this);

        const polyPatch& pp = td.mesh.boundaryMesh()[patchi];

        cloud.functions().postPatch(p, pp);
    }

    return result;
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::MomentumParcel<ParcelType>::hitWallPatch
(
    TrackCloudType&,
    trackingData& td
)
{
    const polyPatch& pp = td.mesh.boundaryMesh()[this->patch(td.mesh)];

    FatalErrorInFunction
        << "Particle " << this->origId() << " hit " << pp.type() << " patch "
        << pp.name() << " at " << this->position(td.mesh)
        << " but no interaction model was specified for this patch"
        << exit(FatalError);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::MomentumParcel<ParcelType>::hitBasicPatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::hitBasicPatch(cloud, td);

    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);

    const polyPatch& pp = td.mesh.boundaryMesh()[p.patch(td.mesh)];

    cloud.functions().postPatch(p, pp);
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "MomentumParcelIO.C"

// ************************************************************************* //
