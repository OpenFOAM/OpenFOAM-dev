/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "SprayParcel.H"
#include "CompositionModel.H"
#include "AtomizationModel.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    ParcelType::setCellValues(td, dt, cellI);
}


template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::cellValueSourceCorrection
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    ParcelType::cellValueSourceCorrection(td, dt, cellI);
}


template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::calc
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    typedef typename TrackData::cloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        td.cloud().composition();

    // Check if parcel belongs to liquid core
    if (liquidCore() > 0.5)
    {
        // Liquid core parcels should not experience coupled forces
        td.cloud().forces().setCalcCoupled(false);
    }

    // Get old mixture composition
    const scalarField& Y0(this->Y());
    scalarField X0(composition.liquids().X(Y0));

    // Check if we have critical or boiling conditions
    scalar TMax = composition.liquids().Tc(X0);
    const scalar T0 = this->T();
    const scalar pc0 = this->pc_;
    if (composition.liquids().pv(pc0, T0, X0) >= pc0*0.999)
    {
        // Set TMax to boiling temperature
        TMax = composition.liquids().pvInvert(pc0, X0);
    }

    // Set the maximum temperature limit
    td.cloud().constProps().setTMax(TMax);

    // Store the parcel properties
    this->Cp() = composition.liquids().Cp(pc0, T0, X0);
    sigma_ = composition.liquids().sigma(pc0, T0, X0);
    const scalar rho0 = composition.liquids().rho(pc0, T0, X0);
    this->rho() = rho0;
    const scalar mass0 = this->mass();
    mu_ = composition.liquids().mu(pc0, T0, X0);

    ParcelType::calc(td, dt, cellI);

    if (td.keepParticle)
    {
        // Reduce the stripped parcel mass due to evaporation
        // assuming the number of particles remains unchanged
        this->ms() -= this->ms()*(mass0 - this->mass())/mass0;

        // Update Cp, sigma, density and diameter due to change in temperature
        // and/or composition
        scalar T1 = this->T();
        const scalarField& Y1(this->Y());
        scalarField X1(composition.liquids().X(Y1));

        this->Cp() = composition.liquids().Cp(this->pc_, T1, X1);

        sigma_ = composition.liquids().sigma(this->pc_, T1, X1);

        scalar rho1 = composition.liquids().rho(this->pc_, T1, X1);
        this->rho() = rho1;

        mu_ = composition.liquids().mu(this->pc_, T1, X1);

        scalar d1 = this->d()*cbrt(rho0/rho1);
        this->d() = d1;

        if (liquidCore() > 0.5)
        {
            calcAtomization(td, dt, cellI);

            // Preserve the total mass/volume by increasing the number of
            // particles in parcels due to breakup
            scalar d2 = this->d();
            this->nParticle() *= pow3(d1/d2);
        }
        else
        {
            calcBreakup(td, dt, cellI);
        }
    }

    // Restore coupled forces
    td.cloud().forces().setCalcCoupled(true);
}


template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::calcAtomization
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    typedef typename TrackData::cloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        td.cloud().composition();

    typedef typename TrackData::cloudType::sprayCloudType sprayCloudType;
    const AtomizationModel<sprayCloudType>& atomization =
        td.cloud().atomization();

    // Average molecular weight of carrier mix - assumes perfect gas
    scalar Wc = this->rhoc_*specie::RR*this->Tc()/this->pc();
    scalar R = specie::RR/Wc;
    scalar Tav = atomization.Taverage(this->T(), this->Tc());

    // Calculate average gas density based on average temperature
    scalar rhoAv = this->pc()/(R*Tav);

    scalar soi = td.cloud().injectors().timeStart();
    scalar currentTime = td.cloud().db().time().value();
    const vector& pos = this->position();
    const vector& injectionPos = this->position0();

    // Disregard the continous phase when calculating the relative velocity
    // (in line with the deactivated coupled assumption)
    scalar Urel = mag(this->U());

    scalar t0 = max(0.0, currentTime - this->age() - soi);
    scalar t1 = min(t0 + dt, td.cloud().injectors().timeEnd() - soi);

    // This should be the vol flow rate from when the parcel was injected
    scalar volFlowRate = td.cloud().injectors().volumeToInject(t0, t1)/dt;

    scalar chi = 0.0;
    if (atomization.calcChi())
    {
        chi = this->chi(td, composition.liquids().X(this->Y()));
    }

    atomization.update
    (
        dt,
        this->d(),
        this->liquidCore(),
        this->tc(),
        this->rho(),
        mu_,
        sigma_,
        volFlowRate,
        rhoAv,
        Urel,
        pos,
        injectionPos,
        td.cloud().pAmbient(),
        chi,
        td.cloud().rndGen()
    );
}


template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::calcBreakup
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    typedef typename TrackData::cloudType cloudType;
    typedef typename cloudType::parcelType parcelType;
    typedef typename cloudType::forceType forceType;

    const parcelType& p = static_cast<const parcelType&>(*this);
    const forceType& forces = td.cloud().forces();

    if (td.cloud().breakup().solveOscillationEq())
    {
        solveTABEq(td, dt);
    }

    // Average molecular weight of carrier mix - assumes perfect gas
    scalar Wc = this->rhoc()*specie::RR*this->Tc()/this->pc();
    scalar R = specie::RR/Wc;
    scalar Tav = td.cloud().atomization().Taverage(this->T(), this->Tc());

    // Calculate average gas density based on average temperature
    scalar rhoAv = this->pc()/(R*Tav);
    scalar muAv = this->muc();
    vector Urel = this->U() - this->Uc();
    scalar Urmag = mag(Urel);
    scalar Re = this->Re(this->U(), this->d(), rhoAv, muAv);

    const scalar mass = p.mass();
    const forceSuSp Fcp = forces.calcCoupled(p, dt, mass, Re, muAv);
    const forceSuSp Fncp = forces.calcNonCoupled(p, dt, mass, Re, muAv);
    this->tMom() = mass/(Fcp.Sp() + Fncp.Sp());

    const vector g = td.cloud().g().value();

    scalar parcelMassChild = 0.0;
    scalar dChild = 0.0;
    if
    (
        td.cloud().breakup().update
        (
            dt,
            g,
            this->d(),
            this->tc(),
            this->ms(),
            this->nParticle(),
            this->KHindex(),
            this->y(),
            this->yDot(),
            this->d0(),
            this->rho(),
            mu_,
            sigma_,
            this->U(),
            rhoAv,
            muAv,
            Urel,
            Urmag,
            this->tMom(),
            dChild,
            parcelMassChild
        )
    )
    {
        scalar Re = rhoAv*Urmag*dChild/muAv;

        // Add child parcel as copy of parent
        SprayParcel<ParcelType>* child = new SprayParcel<ParcelType>(*this);
        child->d() = dChild;
        child->d0() = dChild;
        const scalar massChild = child->mass();
        child->mass0() = massChild;
        child->nParticle() = parcelMassChild/massChild;

        const forceSuSp Fcp =
            forces.calcCoupled(*child, dt, massChild, Re, muAv);
        const forceSuSp Fncp =
            forces.calcNonCoupled(*child, dt, massChild, Re, muAv);

        child->age() = 0.0;
        child->liquidCore() = 0.0;
        child->KHindex() = 1.0;
        child->y() = td.cloud().breakup().y0();
        child->yDot() = td.cloud().breakup().yDot0();
        child->tc() = 0.0;
        child->ms() = -GREAT;
        child->injector() = this->injector();
        child->tMom() = massChild/(Fcp.Sp() + Fncp.Sp());
        child->user() = 0.0;
        child->setCellValues(td, dt, cellI);

        td.cloud().addParticle(child);
    }
}


template<class ParcelType>
template<class TrackData>
Foam::scalar Foam::SprayParcel<ParcelType>::chi
(
    TrackData& td,
    const scalarField& X
) const
{
    // Modifications to take account of the flash boiling on primary break-up

    typedef typename TrackData::cloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        td.cloud().composition();

    scalar chi = 0.0;
    scalar T0 = this->T();
    scalar p0 = this->pc();
    scalar pAmb = td.cloud().pAmbient();

    scalar pv = composition.liquids().pv(p0, T0, X);

    forAll(composition.liquids(), i)
    {
        if (pv >= 0.999*pAmb)
        {
            // Liquid is boiling - calc boiling temperature

            const liquidProperties& liq = composition.liquids().properties()[i];
            scalar TBoil = liq.pvInvert(p0);

            scalar hl = liq.hl(pAmb, TBoil);
            scalar iTp = liq.h(pAmb, T0) - pAmb/liq.rho(pAmb, T0);
            scalar iTb = liq.h(pAmb, TBoil) - pAmb/liq.rho(pAmb, TBoil);

            chi += X[i]*(iTp - iTb)/hl;
        }
    }

    chi = min(1.0, max(chi, 0.0));

    return chi;
}


template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::solveTABEq
(
    TrackData& td,
    const scalar dt
)
{
    const scalar& TABCmu = td.cloud().breakup().TABCmu();
    const scalar& TABWeCrit = td.cloud().breakup().TABWeCrit();
    const scalar& TABComega = td.cloud().breakup().TABComega();

    scalar r = 0.5*this->d();
    scalar r2 = r*r;
    scalar r3 = r*r2;

    // Inverse of characteristic viscous damping time
    scalar rtd = 0.5*TABCmu*mu_/(this->rho()*r2);

    // Oscillation frequency (squared)
    scalar omega2 = TABComega*sigma_/(this->rho()*r3) - rtd*rtd;

    if (omega2 > 0)
    {
        scalar omega = sqrt(omega2);
        scalar rhoc = this->rhoc();
        scalar We = this->We(this->U(), r, rhoc, sigma_)/TABWeCrit;

        scalar y1 = this->y() - We;
        scalar y2 = this->yDot()/omega;

        // Update distortion parameters
        scalar c = cos(omega*dt);
        scalar s = sin(omega*dt);
        scalar e = exp(-rtd*dt);
        y2 = (this->yDot() + y1*rtd)/omega;

        this->y() = We + e*(y1*c + y2*s);
        if (this->y() < 0)
        {
            this->y() = 0.0;
            this->yDot() = 0.0;
        }
        else
        {
            this->yDot() = (We - this->y())*rtd + e*omega*(y2*c - y1*s);
        }
    }
    else
    {
        // Reset distortion parameters
        this->y() = 0;
        this->yDot() = 0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::SprayParcel<ParcelType>::SprayParcel(const SprayParcel<ParcelType>& p)
:
    ParcelType(p),
    d0_(p.d0_),
    position0_(p.position0_),
    sigma_(p.sigma_),
    mu_(p.mu_),
    liquidCore_(p.liquidCore_),
    KHindex_(p.KHindex_),
    y_(p.y_),
    yDot_(p.yDot_),
    tc_(p.tc_),
    ms_(p.ms_),
    injector_(p.injector_),
    tMom_(p.tMom_),
    user_(p.user_)
{}


template<class ParcelType>
Foam::SprayParcel<ParcelType>::SprayParcel
(
    const SprayParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    d0_(p.d0_),
    position0_(p.position0_),
    sigma_(p.sigma_),
    mu_(p.mu_),
    liquidCore_(p.liquidCore_),
    KHindex_(p.KHindex_),
    y_(p.y_),
    yDot_(p.yDot_),
    tc_(p.tc_),
    ms_(p.ms_),
    injector_(p.injector_),
    tMom_(p.tMom_),
    user_(p.user_)
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "SprayParcelIO.C"


// ************************************************************************* //
