/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "StandardChemistryModel.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::StandardChemistryModel<ThermoType>::StandardChemistryModel
(
    const fluidReactionThermo& thermo
)
:
    basicChemistryModel(thermo),
    ODESystem(),
    Y_(this->thermo().composition().Y()),
    mixture_(refCast<const multiComponentMixture<ThermoType>>(this->thermo())),
    specieThermos_(mixture_.specieThermos()),
    reactions_(mixture_.species(), specieThermos_, this->mesh(), *this),
    nSpecie_(Y_.size()),
    nReaction_(reactions_.size()),
    Treact_(basicChemistryModel::template lookupOrDefault<scalar>("Treact", 0)),
    RR_(nSpecie_),
    c_(nSpecie_),
    dcdt_(nSpecie_)
{
    // Create the fields for the chemistry sources
    forAll(RR_, fieldi)
    {
        RR_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RR." + Y_[fieldi].name(),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermo.T().mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0)
            )
        );
    }

    Info<< "StandardChemistryModel: Number of species = " << nSpecie_
        << " and reactions = " << nReaction_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::StandardChemistryModel<ThermoType>::~StandardChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::StandardChemistryModel<ThermoType>::omega
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcdt
) const
{
    dcdt = Zero;

    forAll(reactions_, i)
    {
        const Reaction<ThermoType>& R = reactions_[i];

        R.omega(p, T, c, li, dcdt);
    }
}


template<class ThermoType>
void Foam::StandardChemistryModel<ThermoType>::derivatives
(
    const scalar t,
    const scalarField& c,
    const label li,
    scalarField& dcdt
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    forAll(c_, i)
    {
        c_[i] = max(c[i], 0);
    }

    dcdt = Zero;

    // Evaluate contributions from reactions
    omega(p, T, c_, li, dcdt);

    // Evaluate the effect on the thermodynamic system ...

    // c*Cp
    scalar ccp = 0;
    for (label i = 0; i < nSpecie_; i++)
    {
        ccp += c_[i]*specieThermos_[i].cp(p, T);
    }

    // dT/dt
    scalar& dTdt = dcdt[nSpecie_];
    for (label i = 0; i < nSpecie_; i++)
    {
        dTdt -= dcdt[i]*specieThermos_[i].ha(p, T);
    }
    dTdt /= ccp;

    // dp/dt = 0 (pressure is assumed constant)
}


template<class ThermoType>
void Foam::StandardChemistryModel<ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    const label li,
    scalarField& dcdt,
    scalarSquareMatrix& J
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    forAll(c_, i)
    {
        c_[i] = max(c[i], 0);
    }

    dcdt = Zero;
    J = Zero;

    // Evaluate contributions from reactions
    forAll(reactions_, ri)
    {
        const Reaction<ThermoType>& R = reactions_[ri];
        scalar omegaI, kfwd, kbwd;
        const labelList null;
        R.dwdc(p, T, c_, li, J, dcdt, omegaI, kfwd, kbwd, false, null);
        R.dwdT(p, T, c_, li, omegaI, kfwd, kbwd, J, false, null, nSpecie_);
    }

    // Evaluate the effect on the thermodynamic system ...

    // c*Cp
    scalar ccp = 0, dccpdT = 0;
    for (label i = 0; i < nSpecie_; i++)
    {
        ccp += c_[i]*specieThermos_[i].cp(p, T);
        dccpdT += c_[i]*specieThermos_[i].dcpdT(p, T);
    }

    // dT/dt
    scalar& dTdt = dcdt[nSpecie_];
    for (label i = 0; i < nSpecie_; i++)
    {
        dTdt -= dcdt[i]*specieThermos_[i].ha(p, T);
    }
    dTdt /= ccp;

    // dp/dt = 0 (pressure is assumed constant)

    // d(dTdt)/dc
    for (label i = 0; i < nSpecie_; i++)
    {
        scalar& d2Tdtdci = J(nSpecie_, i);
        for (label j = 0; j < nSpecie_; j++)
        {
            const scalar d2cjdtdci = J(j, i);
            d2Tdtdci -= d2cjdtdci*specieThermos_[j].ha(p, T);
        }
        d2Tdtdci -= specieThermos_[i].cp(p, T)*dTdt;
        d2Tdtdci /= ccp;
    }

    // d(dTdt)/dT
    scalar& d2TdtdT = J(nSpecie_, nSpecie_);
    for (label i = 0; i < nSpecie_; i++)
    {
        const scalar d2cidtdT = J(i, nSpecie_);
        d2TdtdT -=
            dcdt[i]*specieThermos_[i].cp(p, T)
          + d2cidtdT*specieThermos_[i].ha(p, T);
    }
    d2TdtdT -= dTdt*dccpdT;
    d2TdtdT /= ccp;

    // d(dpdt)/dc = 0 (pressure is assumed constant)

    // d(dpdt)/dT = 0 (pressure is assumed constant)
}


template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::StandardChemistryModel<ThermoType>::tc() const
{
    tmp<volScalarField> ttc
    (
        volScalarField::New
        (
            "tc",
            this->mesh(),
            dimensionedScalar(dimTime, small),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    scalarField& tc = ttc.ref();

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    const label nReaction = reactions_.size();

    scalar omegaf, omegar;

    if (this->chemistry_)
    {
        reactionEvaluationScope scope(*this);

        forAll(rho, celli)
        {
            const scalar rhoi = rho[celli];
            const scalar Ti = T[celli];
            const scalar pi = p[celli];

            scalar cSum = 0;

            for (label i=0; i<nSpecie_; i++)
            {
                c_[i] = rhoi*Y_[i][celli]/specieThermos_[i].W();
                cSum += c_[i];
            }

            forAll(reactions_, i)
            {
                const Reaction<ThermoType>& R = reactions_[i];
                R.omega(pi, Ti, c_, celli, omegaf, omegar);

                forAll(R.rhs(), s)
                {
                    tc[celli] += R.rhs()[s].stoichCoeff*omegaf;
                }
            }

            tc[celli] = nReaction*cSum/tc[celli];
        }
    }

    ttc.ref().correctBoundaryConditions();

    return ttc;
}


template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::StandardChemistryModel<ThermoType>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        volScalarField::New
        (
            "Qdot",
            this->mesh_,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (this->chemistry_)
    {
        reactionEvaluationScope scope(*this);

        scalarField& Qdot = tQdot.ref();

        forAll(Y_, i)
        {
            forAll(Qdot, celli)
            {
                const scalar hi = specieThermos_[i].Hf();
                Qdot[celli] -= hi*RR_[i][celli];
            }
        }
    }

    return tQdot;
}


template<class ThermoType>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::StandardChemistryModel<ThermoType>::calculateRR
(
    const label ri,
    const label si
) const
{
    tmp<volScalarField::Internal> tRR
    (
        volScalarField::Internal::New
        (
            "RR",
            this->mesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, 0)
        )
    );

    volScalarField::Internal& RR = tRR.ref();

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    reactionEvaluationScope scope(*this);

    scalar omegaf, omegar;

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T[celli];
        const scalar pi = p[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Y_[i][celli];
            c_[i] = rhoi*Yi/specieThermos_[i].W();
        }

        const Reaction<ThermoType>& R = reactions_[ri];
        const scalar omegaI = R.omega(pi, Ti, c_, celli, omegaf, omegar);

        forAll(R.lhs(), s)
        {
            if (si == R.lhs()[s].index)
            {
                RR[celli] -= R.lhs()[s].stoichCoeff*omegaI;
            }
        }

        forAll(R.rhs(), s)
        {
            if (si == R.rhs()[s].index)
            {
                RR[celli] += R.rhs()[s].stoichCoeff*omegaI;
            }
        }

        RR[celli] *= specieThermos_[si].W();
    }

    return tRR;
}


template<class ThermoType>
void Foam::StandardChemistryModel<ThermoType>::calculate()
{
    if (!this->chemistry_)
    {
        return;
    }

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    reactionEvaluationScope scope(*this);

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T[celli];
        const scalar pi = p[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Y_[i][celli];
            c_[i] = rhoi*Yi/specieThermos_[i].W();
        }

        omega(pi, Ti, c_, celli, dcdt_);

        for (label i=0; i<nSpecie_; i++)
        {
            RR_[i][celli] = dcdt_[i]*specieThermos_[i].W();
        }
    }
}


template<class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::StandardChemistryModel<ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    basicChemistryModel::correct();

    scalar deltaTMin = great;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    tmp<volScalarField> trho(this->thermo().rho0());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T().oldTime();
    const scalarField& p = this->thermo().p().oldTime();

    reactionEvaluationScope scope(*this);

    scalarField c0(nSpecie_);

    forAll(rho, celli)
    {
        scalar Ti = T[celli];

        if (Ti > Treact_)
        {
            const scalar rhoi = rho[celli];
            scalar pi = p[celli];

            for (label i=0; i<nSpecie_; i++)
            {
                c_[i] = rhoi*Y_[i].oldTime()[celli]/specieThermos_[i].W();
                c0[i] = c_[i];
            }

            // Initialise time progress
            scalar timeLeft = deltaT[celli];

            // Calculate the chemical source terms
            while (timeLeft > small)
            {
                scalar dt = timeLeft;
                this->solve(pi, Ti, c_, celli, dt, this->deltaTChem_[celli]);
                timeLeft -= dt;
            }

            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

            this->deltaTChem_[celli] =
                min(this->deltaTChem_[celli], this->deltaTChemMax_);

            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] =
                    (c_[i] - c0[i])*specieThermos_[i].W()/deltaT[celli];
            }
        }
        else
        {
            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] = 0;
            }
        }
    }

    return deltaTMin;
}


template<class ThermoType>
Foam::scalar Foam::StandardChemistryModel<ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}


template<class ThermoType>
Foam::scalar Foam::StandardChemistryModel<ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


// ************************************************************************* //
