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

#include "ReactingParcel.H"
#include "specie.H"
#include "CompositionModel.H"
#include "PhaseChangeModel.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::calcPhaseChange
(
    TrackData& td,
    const scalar dt,
    const label cellI,
    const scalar Re,
    const scalar Pr,
    const scalar Ts,
    const scalar nus,
    const scalar d,
    const scalar T,
    const scalar mass,
    const label idPhase,
    const scalar YPhase,
    const scalarField& YComponents,
    scalarField& dMassPC,
    scalar& Sh,
    scalar& N,
    scalar& NCpW,
    scalarField& Cs
)
{
    typedef typename TrackData::cloudType::reactingCloudType reactingCloudType;
    PhaseChangeModel<reactingCloudType>& phaseChange = td.cloud().phaseChange();

    if (!phaseChange.active() || (YPhase < SMALL))
    {
        return;
    }

    scalar Tvap = phaseChange.Tvap(YComponents);

    if (T < Tvap)
    {
        return;
    }

    const scalar TMax = phaseChange.TMax(pc_, YComponents);
    const scalar Tdash = min(T, TMax);
    const scalar Tsdash = min(Ts, TMax);

    // Calculate mass transfer due to phase change
    phaseChange.calculate
    (
        dt,
        cellI,
        Re,
        Pr,
        d,
        nus,
        Tdash,
        Tsdash,
        pc_,
        this->Tc_,
        YComponents,
        dMassPC
    );

    // Limit phase change mass by availability of each specie
    dMassPC = min(mass*YPhase*YComponents, dMassPC);

    const scalar dMassTot = sum(dMassPC);

    // Add to cumulative phase change mass
    phaseChange.addToPhaseChangeMass(this->nParticle_*dMassTot);

    const CompositionModel<reactingCloudType>& composition =
        td.cloud().composition();

    forAll(dMassPC, i)
    {
        const label idc = composition.localToGlobalCarrierId(idPhase, i);
        const label idl = composition.globalIds(idPhase)[i];

        const scalar dh = phaseChange.dh(idc, idl, pc_, Tdash);
        Sh -= dMassPC[i]*dh/dt;
    }


    // Update molar emissions
    if (td.cloud().heatTransfer().BirdCorrection())
    {
        // Average molecular weight of carrier mix - assumes perfect gas
        const scalar Wc = this->rhoc_*specie::RR*this->Tc_/this->pc_;

        forAll(dMassPC, i)
        {
            const label idc = composition.localToGlobalCarrierId(idPhase, i);
            const label idl = composition.globalIds(idPhase)[i];

            const scalar Cp = composition.carrier().Cp(idc, pc_, Tsdash);
            const scalar W = composition.carrier().W(idc);
            const scalar Ni = dMassPC[i]/(this->areaS(d)*dt*W);

            const scalar Dab =
                composition.liquids().properties()[idl].D(pc_, Tsdash, Wc);

            // Molar flux of species coming from the particle (kmol/m^2/s)
            N += Ni;

            // Sum of Ni*Cpi*Wi of emission species
            NCpW += Ni*Cp*W;

            // Concentrations of emission species
            Cs[idc] += Ni*d/(2.0*Dab);
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
    if (mass1 > ROOTVSMALL)
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
    Y_(p.Y_),
    pc_(p.pc_)
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
    Y_(p.Y_),
    pc_(p.pc_)
{}


// * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    ParcelType::setCellValues(td, dt, cellI);

    pc_ = td.pInterp().interpolate
    (
        this->position(),
        this->currentTetIndices()
    );

    if (pc_ < td.cloud().constProps().pMin())
    {
        if (debug)
        {
            WarningIn
            (
                "void Foam::ReactingParcel<ParcelType>::setCellValues"
                "("
                    "TrackData&, "
                    "const scalar, "
                    "const label"
                ")"
            )   << "Limiting observed pressure in cell " << cellI << " to "
                << td.cloud().constProps().pMin() <<  nl << endl;
        }

        pc_ = td.cloud().constProps().pMin();
    }
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::cellValueSourceCorrection
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    scalar addedMass = 0.0;
    scalar maxMassI = 0.0;
    forAll(td.cloud().rhoTrans(), i)
    {
        scalar dm = td.cloud().rhoTrans(i)[cellI];
        maxMassI = max(maxMassI, mag(dm));
        addedMass += dm;
    }

    if (maxMassI < ROOTVSMALL)
    {
        return;
    }

    const scalar massCell = this->massCell(cellI);

    this->rhoc_ += addedMass/td.cloud().pMesh().cellVolumes()[cellI];

    const scalar massCellNew = massCell + addedMass;
    this->Uc_ = (this->Uc_*massCell + td.cloud().UTrans()[cellI])/massCellNew;

    scalar CpEff = 0.0;
    forAll(td.cloud().rhoTrans(), i)
    {
        scalar Y = td.cloud().rhoTrans(i)[cellI]/addedMass;
        CpEff += Y*td.cloud().composition().carrier().Cp
        (
            i,
            this->pc_,
            this->Tc_
        );
    }

    const scalar Cpc = td.CpInterp().psi()[cellI];
    this->Cpc_ = (massCell*Cpc + addedMass*CpEff)/massCellNew;

    this->Tc_ += td.cloud().hsTrans()[cellI]/(this->Cpc_*massCellNew);

    if (this->Tc_ < td.cloud().constProps().TMin())
    {
        if (debug)
        {
            WarningIn
            (
                "void Foam::ReactingParcel<ParcelType>::"
                "cellValueSourceCorrection"
                "("
                    "TrackData&, "
                    "const scalar, "
                    "const label"
                ")"
            )   << "Limiting observed temperature in cell " << cellI << " to "
                << td.cloud().constProps().TMin() <<  nl << endl;
        }

        this->Tc_ = td.cloud().constProps().TMin();
    }

//  constant pressure
//  this->pc_ = this->pc_;
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::correctSurfaceValues
(
    TrackData& td,
    const label cellI,
    const scalar T,
    const scalarField& Cs,
    scalar& rhos,
    scalar& mus,
    scalar& Prs,
    scalar& kappas
)
{
    // No correction if total concentration of emitted species is small
    if (!td.cloud().heatTransfer().BirdCorrection() || (sum(Cs) < SMALL))
    {
        return;
    }

    const SLGThermo& thermo = td.cloud().thermo();

    // Far field carrier  molar fractions
    scalarField Xinf(thermo.carrier().species().size());

    forAll(Xinf, i)
    {
        Xinf[i] = thermo.carrier().Y(i)[cellI]/thermo.carrier().W(i);
    }
    Xinf /= sum(Xinf);

    // Molar fraction of far field species at particle surface
    const scalar Xsff = 1.0 - min(sum(Cs)*specie::RR*this->T_/pc_, 1.0);

    // Surface carrier total molar concentration
    const scalar CsTot = pc_/(specie::RR*this->T_);

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
        mus += Ys[i]*sqrtW*thermo.carrier().mu(i, pc_, T);
        kappas += Ys[i]*cbrtW*thermo.carrier().kappa(i, pc_, T);
        Cps += Xs[i]*thermo.carrier().Cp(i, pc_, T);

        sumYiSqrtW += Ys[i]*sqrtW;
        sumYiCbrtW += Ys[i]*cbrtW;
    }

    Cps = max(Cps, ROOTVSMALL);

    rhos *= pc_/(specie::RR*T);
    rhos = max(rhos, ROOTVSMALL);

    mus /= sumYiSqrtW;
    mus = max(mus, ROOTVSMALL);

    kappas /= sumYiCbrtW;
    kappas = max(kappas, ROOTVSMALL);

    Prs = Cps*mus/kappas;
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::calc
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    typedef typename TrackData::cloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        td.cloud().composition();


    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const scalar np0 = this->nParticle_;
    const scalar d0 = this->d_;
    const vector& U0 = this->U_;
    const scalar T0 = this->T_;
    const scalar mass0 = this->mass();


    // Calc surface values
    scalar Ts, rhos, mus, Prs, kappas;
    this->calcSurfaceValues(td, cellI, T0, Ts, rhos, mus, Prs, kappas);
    scalar Res = this->Re(U0, d0, rhos, mus);


    // Sources
    //~~~~~~~~

    // Explicit momentum source for particle
    vector Su = vector::zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = vector::zero;

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
        td,
        dt,
        cellI,
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

    this->Cp_ = composition.Cp(0, Y_, pc_, T0);

    // Update particle density or diameter
    if (td.cloud().constProps().constantVolume())
    {
        this->rho_ = mass1/this->volume();
    }
    else
    {
        this->d_ = cbrt(mass1/this->rho_*6.0/pi);
    }

    // Remove the particle when mass falls below minimum threshold
    if (np0*mass1 < td.cloud().constProps().minParticleMass())
    {
        td.keepParticle = false;

        if (td.cloud().solution().coupled())
        {
            scalar dm = np0*mass0;

            // Absorb parcel into carrier phase
            forAll(Y_, i)
            {
                scalar dmi = dm*Y_[i];
                label gid = composition.localToGlobalCarrierId(0, i);
                scalar hs = composition.carrier().Hs(gid, pc_, T0);

                td.cloud().rhoTrans(gid)[cellI] += dmi;
                td.cloud().hsTrans()[cellI] += dmi*hs;
            }
            td.cloud().UTrans()[cellI] += dm*U0;

            td.cloud().phaseChange().addToPhaseChangeMass(np0*mass1);
        }

        return;
    }

    // Correct surface values due to emitted species
    correctSurfaceValues(td, cellI, Ts, Cs, rhos, mus, Prs, kappas);
    Res = this->Re(U0, this->d_, rhos, mus);


    // 3. Compute heat- and momentum transfers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Heat transfer
    // ~~~~~~~~~~~~~

    // Calculate new particle temperature
    this->T_ =
        this->calcHeatTransfer
        (
            td,
            dt,
            cellI,
            Res,
            Prs,
            kappas,
            NCpW,
            Sh,
            dhsTrans,
            Sph
        );

    this->Cp_ = composition.Cp(0, Y_, pc_, T0);


    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    this->U_ =
        this->calcVelocity(td, dt, cellI, Res, mus, mass1, Su, dUTrans, Spu);


    // 4. Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (td.cloud().solution().coupled())
    {
        // Transfer mass lost to carrier mass, momentum and enthalpy sources
        forAll(dMass, i)
        {
            scalar dm = np0*dMass[i];
            label gid = composition.localToGlobalCarrierId(0, i);
            scalar hs = composition.carrier().Hs(gid, pc_, T0);

            td.cloud().rhoTrans(gid)[cellI] += dm;
            td.cloud().UTrans()[cellI] += dm*U0;
            td.cloud().hsTrans()[cellI] += dm*hs;
        }

        // Update momentum transfer
        td.cloud().UTrans()[cellI] += np0*dUTrans;
        td.cloud().UCoeff()[cellI] += np0*Spu;

        // Update sensible enthalpy transfer
        td.cloud().hsTrans()[cellI] += np0*dhsTrans;
        td.cloud().hsCoeff()[cellI] += np0*Sph;

        // Update radiation fields
        if (td.cloud().radiation())
        {
            const scalar ap = this->areaP();
            const scalar T4 = pow4(this->T_);
            td.cloud().radAreaP()[cellI] += dt*np0*ap;
            td.cloud().radT4()[cellI] += dt*np0*T4;
            td.cloud().radAreaPT4()[cellI] += dt*np0*ap*T4;
        }
    }
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ReactingParcelIO.C"

// ************************************************************************* //

