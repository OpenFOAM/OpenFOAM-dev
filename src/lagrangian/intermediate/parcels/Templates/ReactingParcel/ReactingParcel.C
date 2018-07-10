/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "ReactingParcel.H"
#include "specie.H"
#include "CompositionModel.H"
#include "PhaseChangeModel.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingParcel<ParcelType>::calcPhaseChange
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar Re,
    const scalar Pr,
    const scalar Ts,
    const scalar nus,
    const scalar d,
    const scalar T,
    const scalar mass,
    const label idPhase,
    const scalar YPhase,
    const scalarField& Y,
    scalarField& dMassPC,
    scalar& Sh,
    scalar& N,
    scalar& NCpW,
    scalarField& Cs
)
{
    typedef typename TrackCloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        cloud.composition();
    PhaseChangeModel<reactingCloudType>& phaseChange = cloud.phaseChange();

    if (!phaseChange.active() || (YPhase < small))
    {
        return;
    }

    scalarField X(composition.liquids().X(Y));

    scalar Tvap = phaseChange.Tvap(X);

    if (T < Tvap)
    {
        return;
    }

    const scalar TMax = phaseChange.TMax(td.pc(), X);
    const scalar Tdash = min(T, TMax);
    const scalar Tsdash = min(Ts, TMax);

    scalarField hmm(dMassPC);

    // Calculate mass transfer due to phase change
    phaseChange.calculate
    (
        dt,
        this->cell(),
        Re,
        Pr,
        d,
        nus,
        Tdash,
        Tsdash,
        td.pc(),
        td.Tc(),
        X,
        dMassPC
    );

    // Limit phase change mass by availability of each specie
    dMassPC = min(mass*YPhase*Y, dMassPC);

    const scalar dMassTot = sum(dMassPC);

    // Add to cumulative phase change mass
    phaseChange.addToPhaseChangeMass(this->nParticle_*dMassTot);

    forAll(dMassPC, i)
    {
        const label cid = composition.localToCarrierId(idPhase, i);

        const scalar dh = phaseChange.dh(cid, i, td.pc(), Tdash);
        Sh -= dMassPC[i]*dh/dt;
    }


    // Update molar emissions
    if (cloud.heatTransfer().BirdCorrection())
    {
        // Average molecular weight of carrier mix - assumes perfect gas
        const scalar Wc = td.rhoc()*RR*td.Tc()/td.pc();

        forAll(dMassPC, i)
        {
            const label cid = composition.localToCarrierId(idPhase, i);

            const scalar Cp = composition.carrier().Cp(cid, td.pc(), Tsdash);
            const scalar W = composition.carrier().W(cid);
            const scalar Ni = dMassPC[i]/(this->areaS(d)*dt*W);

            const scalar Dab =
                composition.liquids().properties()[i].D(td.pc(), Tsdash, Wc);

            // Molar flux of species coming from the particle (kmol/m^2/s)
            N += Ni;

            // Sum of Ni*Cpi*Wi of emission species
            NCpW += Ni*Cp*W;

            // Concentrations of emission species
            Cs[cid] += Ni*d/(2.0*Dab);
        }
    }
}


template<class ParcelType>
Foam::scalar Foam::ReactingParcel<ParcelType>::updateMassFraction
(
    const scalar mass0,
    const scalarField& dMass,
    scalarField& Y
) const
{
    scalar mass1 = mass0 - sum(dMass);

    // only update the mass fractions if the new particle mass is finite
    if (mass1 > rootVSmall)
    {
        forAll(Y, i)
        {
            Y[i] = (Y[i]*mass0 - dMass[i])/mass1;
        }
    }

    return mass1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingParcel<ParcelType>::ReactingParcel
(
    const ReactingParcel<ParcelType>& p
)
:
    ParcelType(p),
    mass0_(p.mass0_),
    Y_(p.Y_)
{}


template<class ParcelType>
Foam::ReactingParcel<ParcelType>::ReactingParcel
(
    const ReactingParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    mass0_(p.mass0_),
    Y_(p.Y_)
{}


// * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);

    td.pc() = td.pInterp().interpolate
    (
        this->coordinates(),
        this->currentTetIndices()
    );

    if (td.pc() < cloud.constProps().pMin())
    {
        if (debug)
        {
            WarningInFunction
                << "Limiting observed pressure in cell " << this->cell()
                << " to " << cloud.constProps().pMin() <<  nl << endl;
        }

        td.pc() = cloud.constProps().pMin();
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    scalar addedMass = 0.0;
    scalar maxMassI = 0.0;
    forAll(cloud.rhoTrans(), i)
    {
        scalar dm = cloud.rhoTrans(i)[this->cell()];
        maxMassI = max(maxMassI, mag(dm));
        addedMass += dm;
    }

    if (maxMassI < rootVSmall)
    {
        return;
    }

    const scalar massCell = this->massCell(td);

    td.rhoc() += addedMass/cloud.pMesh().cellVolumes()[this->cell()];

    const scalar massCellNew = massCell + addedMass;
    td.Uc() = (td.Uc()*massCell + cloud.UTrans()[this->cell()])/massCellNew;

    scalar CpEff = 0.0;
    forAll(cloud.rhoTrans(), i)
    {
        scalar Y = cloud.rhoTrans(i)[this->cell()]/addedMass;
        CpEff += Y*cloud.composition().carrier().Cp(i, td.pc(), td.Tc());
    }

    const scalar Cpc = td.CpInterp().psi()[this->cell()];
    td.Cpc() = (massCell*Cpc + addedMass*CpEff)/massCellNew;

    td.Tc() += cloud.hsTrans()[this->cell()]/(td.Cpc()*massCellNew);

    if (td.Tc() < cloud.constProps().TMin())
    {
        if (debug)
        {
            WarningInFunction
                << "Limiting observed temperature in cell " << this->cell()
                << " to " << cloud.constProps().TMin() <<  nl << endl;
        }

        td.Tc() = cloud.constProps().TMin();
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingParcel<ParcelType>::correctSurfaceValues
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar T,
    const scalarField& Cs,
    scalar& rhos,
    scalar& mus,
    scalar& Prs,
    scalar& kappas
)
{
    // No correction if total concentration of emitted species is small
    if (!cloud.heatTransfer().BirdCorrection() || (sum(Cs) < small))
    {
        return;
    }

    const SLGThermo& thermo = cloud.thermo();

    // Far field carrier  molar fractions
    scalarField Xinf(thermo.carrier().species().size());

    forAll(Xinf, i)
    {
        Xinf[i] = thermo.carrier().Y(i)[this->cell()]/thermo.carrier().W(i);
    }
    Xinf /= sum(Xinf);

    // Molar fraction of far field species at particle surface
    const scalar Xsff = 1.0 - min(sum(Cs)*RR*this->T_/td.pc(), 1.0);

    // Surface carrier total molar concentration
    const scalar CsTot = td.pc()/(RR*this->T_);

    // Surface carrier composition (molar fraction)
    scalarField Xs(Xinf.size());

    // Surface carrier composition (mass fraction)
    scalarField Ys(Xinf.size());

    forAll(Xs, i)
    {
        // Molar concentration of species at particle surface
        const scalar Csi = Cs[i] + Xsff*Xinf[i]*CsTot;

        Xs[i] = (2.0*Csi + Xinf[i]*CsTot)/3.0;
        Ys[i] = Xs[i]*thermo.carrier().W(i);
    }
    Xs /= sum(Xs);
    Ys /= sum(Ys);


    rhos = 0;
    mus = 0;
    kappas = 0;
    scalar Cps = 0;
    scalar sumYiSqrtW = 0;
    scalar sumYiCbrtW = 0;

    forAll(Ys, i)
    {
        const scalar W = thermo.carrier().W(i);
        const scalar sqrtW = sqrt(W);
        const scalar cbrtW = cbrt(W);

        rhos += Xs[i]*W;
        mus += Ys[i]*sqrtW*thermo.carrier().mu(i, td.pc(), T);
        kappas += Ys[i]*cbrtW*thermo.carrier().kappa(i, td.pc(), T);
        Cps += Xs[i]*thermo.carrier().Cp(i, td.pc(), T);

        sumYiSqrtW += Ys[i]*sqrtW;
        sumYiCbrtW += Ys[i]*cbrtW;
    }

    Cps = max(Cps, rootVSmall);

    rhos *= td.pc()/(RR*T);
    rhos = max(rhos, rootVSmall);

    mus /= sumYiSqrtW;
    mus = max(mus, rootVSmall);

    kappas /= sumYiCbrtW;
    kappas = max(kappas, rootVSmall);

    Prs = Cps*mus/kappas;
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    typedef typename TrackCloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        cloud.composition();


    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const scalar np0 = this->nParticle_;
    const scalar d0 = this->d_;
    const vector& U0 = this->U_;
    const scalar T0 = this->T_;
    const scalar mass0 = this->mass();


    // Calc surface values
    scalar Ts, rhos, mus, Prs, kappas;
    this->calcSurfaceValues(cloud, td, T0, Ts, rhos, mus, Prs, kappas);
    scalar Res = this->Re(rhos, U0, td.Uc(), d0, mus);


    // Sources
    // ~~~~~~~

    // Explicit momentum source for particle
    vector Su = Zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = Zero;

    // Explicit enthalpy source for particle
    scalar Sh = 0.0;

    // Linearised enthalpy source coefficient
    scalar Sph = 0.0;

    // Sensible enthalpy transfer from the particle to the carrier phase
    scalar dhsTrans = 0.0;


    // 1. Compute models that contribute to mass transfer - U, T held constant
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Phase change
    // ~~~~~~~~~~~~

    // Mass transfer due to phase change
    scalarField dMassPC(Y_.size(), 0.0);

    // Molar flux of species emitted from the particle (kmol/m^2/s)
    scalar Ne = 0.0;

    // Sum Ni*Cpi*Wi of emission species
    scalar NCpW = 0.0;

    // Surface concentrations of emitted species
    scalarField Cs(composition.carrier().species().size(), 0.0);

    // Calc mass and enthalpy transfer due to phase change
    calcPhaseChange
    (
        cloud,
        td,
        dt,
        Res,
        Prs,
        Ts,
        mus/rhos,
        d0,
        T0,
        mass0,
        0,
        1.0,
        Y_,
        dMassPC,
        Sh,
        Ne,
        NCpW,
        Cs
    );


    // 2. Update the parcel properties due to change in mass
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalarField dMass(dMassPC);
    scalar mass1 = updateMassFraction(mass0, dMass, Y_);

    this->Cp_ = composition.Cp(0, Y_, td.pc(), T0);

    // Update particle density or diameter
    if (cloud.constProps().constantVolume())
    {
        this->rho_ = mass1/this->volume();
    }
    else
    {
        this->d_ = cbrt(mass1/this->rho_*6.0/pi);
    }

    // Remove the particle when mass falls below minimum threshold
    if (np0*mass1 < cloud.constProps().minParcelMass())
    {
        td.keepParticle = false;

        if (cloud.solution().coupled())
        {
            scalar dm = np0*mass0;

            // Absorb parcel into carrier phase
            forAll(Y_, i)
            {
                scalar dmi = dm*Y_[i];
                label gid = composition.localToCarrierId(0, i);
                scalar hs = composition.carrier().Hs(gid, td.pc(), T0);

                cloud.rhoTrans(gid)[this->cell()] += dmi;
                cloud.hsTrans()[this->cell()] += dmi*hs;
            }
            cloud.UTrans()[this->cell()] += dm*U0;

            cloud.phaseChange().addToPhaseChangeMass(np0*mass1);
        }

        return;
    }

    // Correct surface values due to emitted species
    correctSurfaceValues(cloud, td, Ts, Cs, rhos, mus, Prs, kappas);
    Res = this->Re(rhos, U0, td.Uc(), this->d(), mus);


    // 3. Compute heat- and momentum transfers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Heat transfer
    // ~~~~~~~~~~~~~

    // Calculate new particle temperature
    this->T_ =
        this->calcHeatTransfer
        (
            cloud,
            td,
            dt,
            Res,
            Prs,
            kappas,
            NCpW,
            Sh,
            dhsTrans,
            Sph
        );

    this->Cp_ = composition.Cp(0, Y_, td.pc(), T0);


    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    this->U_ =
        this->calcVelocity(cloud, td, dt, Res, mus, mass1, Su, dUTrans, Spu);


    // 4. Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (cloud.solution().coupled())
    {
        // Transfer mass lost to carrier mass, momentum and enthalpy sources
        forAll(dMass, i)
        {
            scalar dm = np0*dMass[i];
            label gid = composition.localToCarrierId(0, i);
            scalar hs = composition.carrier().Hs(gid, td.pc(), T0);

            cloud.rhoTrans(gid)[this->cell()] += dm;
            cloud.UTrans()[this->cell()] += dm*U0;
            cloud.hsTrans()[this->cell()] += dm*hs;
        }

        // Update momentum transfer
        cloud.UTrans()[this->cell()] += np0*dUTrans;
        cloud.UCoeff()[this->cell()] += np0*Spu;

        // Update sensible enthalpy transfer
        cloud.hsTrans()[this->cell()] += np0*dhsTrans;
        cloud.hsCoeff()[this->cell()] += np0*Sph;

        // Update radiation fields
        if (cloud.radiation())
        {
            const scalar ap = this->areaP();
            const scalar T4 = pow4(T0);
            cloud.radAreaP()[this->cell()] += dt*np0*ap;
            cloud.radT4()[this->cell()] += dt*np0*T4;
            cloud.radAreaPT4()[this->cell()] += dt*np0*ap*T4;
        }
    }
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ReactingParcelIO.C"

// ************************************************************************* //
