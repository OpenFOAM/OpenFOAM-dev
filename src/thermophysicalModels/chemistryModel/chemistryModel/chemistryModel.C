/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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

#include "chemistryModel.H"
#include "UniformField.H"
#include "localEulerDdtScheme.H"
#include "clockTime.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::chemistryModel<ThermoType>::chemistryModel
(
    const fluidReactionThermo& thermo
)
:
    basicChemistryModel(thermo),
    ODESystem(),
    log_(this->lookupOrDefault("log", false)),
    Y_(this->thermo().composition().Y()),
    mixture_(refCast<const multiComponentMixture<ThermoType>>(this->thermo())),
    specieThermos_(mixture_.specieThermos()),
    reactions_(mixture_.species(), specieThermos_, this->mesh(), *this),
    nSpecie_(Y_.size()),
    nReaction_(reactions_.size()),
    Treact_(basicChemistryModel::template lookupOrDefault<scalar>("Treact", 0)),
    RR_(nSpecie_),
    c_(nSpecie_),
    dcdt_(nSpecie_),
    cTos_(nSpecie_, -1),
    sToc_(nSpecie_),
    mechRedPtr_
    (
        chemistryReductionMethod<ThermoType>::New
        (
            *this,
            *this
        )
    ),
    mechRed_(*mechRedPtr_),
    mechRedActive_(mechRed_.active()),
    tabulationPtr_
    (
        chemistryTabulationMethod<ThermoType>::New
        (
            *this,
            *this
        )
    ),
    tabulation_(*tabulationPtr_)
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

    Info<< "chemistryModel: Number of species = " << nSpecie_
        << " and reactions = " << nReaction_ << endl;

    // When the mechanism reduction method is used, the 'active' flag for every
    // species should be initialised (by default 'active' is true)
    if (mechRedActive_)
    {
        const basicSpecieMixture& composition = this->thermo().composition();

        forAll(Y_, i)
        {
            typeIOobject<volScalarField> header
            (
                Y_[i].name(),
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ
            );

            // Check if the species file is provided, if not set inactive
            // and NO_WRITE
            if (!header.headerOk())
            {
                composition.setInactive(i);
            }
        }
    }

    if (log_)
    {
        cpuSolveFile_ = logFile("cpu_solve.out");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::chemistryModel<ThermoType>::~chemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::chemistryModel<ThermoType>::omega
(
    const scalar p,
    const scalar T,
    const scalarField& c, // Contains all species even when mechRed is active
    const label li,
    scalarField& dcdt
) const
{
    if (mechRedActive_)
    {
        forAll(sToc_, si)
        {
            dcdt_[sToc_[si]] = 0;
        }

        forAll(reactions_, i)
        {
            if (!mechRed_.reactionDisabled(i))
            {
                const Reaction<ThermoType>& R = reactions_[i];

                R.omega(p, T, c, li, dcdt_);
            }
        }

        forAll(sToc_, si)
        {
            dcdt[si] = dcdt_[sToc_[si]];
        }
    }
    else
    {
        dcdt = Zero;

        forAll(reactions_, i)
        {
            const Reaction<ThermoType>& R = reactions_[i];
            R.omega(p, T, c, li, dcdt);
        }
    }
}


template<class ThermoType>
void Foam::chemistryModel<ThermoType>::derivatives
(
    const scalar time,
    const scalarField& c,
    const label li,
    scalarField& dcdt
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    if (mechRedActive_)
    {
        forAll(sToc_, i)
        {
            c_[sToc_[i]] = max(c[i], 0);
        }
    }
    else
    {
        forAll(c_, i)
        {
            c_[i] = max(c[i], 0);
        }
    }

    // Evaluate contributions from reactions
    omega(p, T, c_, li, dcdt);

    // Evaluate the effect on the thermodynamic system ...

    // c*Cp
    scalar ccp = 0;
    for (label i=0; i<c_.size(); i++)
    {
        ccp += c_[i]*specieThermos_[i].cp(p, T);
    }

    // dT/dt
    scalar& dTdt = dcdt[nSpecie_];
    dTdt = 0;
    for (label i=0; i<nSpecie_; i++)
    {
        dTdt -= dcdt[i]*specieThermos_[sToc(i)].ha(p, T);
    }
    dTdt /= ccp;

    // dp/dt = 0 (pressure is assumed constant)
    dcdt[nSpecie_ + 1] = 0;
}


template<class ThermoType>
void Foam::chemistryModel<ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    const label li,
    scalarField& dcdt,
    scalarSquareMatrix& J
) const
{
    // If the mechanism reduction is active, the computed Jacobian
    // is compact (size of the reduced set of species)
    // but according to the information of the complete set
    // (i.e. for the third-body efficiencies)

    if (mechRedActive_)
    {
        forAll(sToc_, i)
        {
            c_[sToc_[i]] = max(c[i], 0);
        }
    }
    else
    {
        forAll(c_, i)
        {
            c_[i] = max(c[i], 0);
        }
    }

    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    dcdt = Zero;
    J = Zero;

    // Evaluate contributions from reactions
    forAll(reactions_, ri)
    {
        if (!mechRed_.reactionDisabled(ri))
        {
            const Reaction<ThermoType>& R = reactions_[ri];
            scalar omegaI, kfwd, kbwd;
            R.dwdc
            (
                p,
                T,
                c_,
                li,
                J,
                dcdt,
                omegaI,
                kfwd,
                kbwd,
                mechRedActive_,
                cTos_
            );
            R.dwdT
            (
                p,
                T,
                c_,
                li,
                omegaI,
                kfwd,
                kbwd,
                J,
                mechRedActive_,
                cTos_,
                nSpecie_
            );
        }
    }

    // Evaluate the effect on the thermodynamic system ...

    // c*Cp
    scalar ccp = 0, dccpdT = 0;
    forAll(c_, i)
    {
        ccp += c_[i]*specieThermos_[i].cp(p, T);
        dccpdT += c_[i]*specieThermos_[i].dcpdT(p, T);
    }

    // dT/dt
    scalar& dTdt = dcdt[nSpecie_];
    for (label i=0; i<nSpecie_; i++)
    {
        dTdt -= dcdt[i]*specieThermos_[sToc(i)].ha(p, T);
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
            d2Tdtdci -= d2cjdtdci*specieThermos_[sToc(j)].ha(p, T);
        }
        d2Tdtdci -= specieThermos_[sToc(i)].cp(p, T)*dTdt;
        d2Tdtdci /= ccp;
    }

    // d(dTdt)/dT
    scalar& d2TdtdT = J(nSpecie_, nSpecie_);
    for (label i = 0; i < nSpecie_; i++)
    {
        const scalar d2cidtdT = J(i, nSpecie_);
        const label si = sToc(i);
        d2TdtdT -=
            dcdt[i]*specieThermos_[si].cp(p, T)
          + d2cidtdT*specieThermos_[si].ha(p, T);
    }
    d2TdtdT -= dTdt*dccpdT;
    d2TdtdT /= ccp;

    // d(dpdt)/dc = 0 (pressure is assumed constant)

    // d(dpdt)/dT = 0 (pressure is assumed constant)
}


template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::chemistryModel<ThermoType>::tc() const
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

    if (this->chemistry_)
    {
        reactionEvaluationScope scope(*this);

        forAll(rho, celli)
        {
            const scalar rhoi = rho[celli];
            const scalar Ti = T[celli];
            const scalar pi = p[celli];

            for (label i=0; i<nSpecie_; i++)
            {
                c_[i] = rhoi*Y_[i][celli]/specieThermos_[i].W();
            }

            // A reaction's rate scale is calculated as it's molar
            // production rate divided by the total number of moles in the
            // system.
            //
            // The system rate scale is the average of the reactions' rate
            // scales weighted by the reactions' molar production rates. This
            // weighting ensures that dominant reactions provide the largest
            // contribution to the system rate scale.
            //
            // The system time scale is then the reciprocal of the system rate
            // scale.
            //
            // Contributions from forward and reverse reaction rates are
            // handled independently and identically so that reversible
            // reactions produce the same result as the equivalent pair of
            // irreversible reactions.

            scalar sumW = 0, sumWRateByCTot = 0;
            forAll(reactions_, i)
            {
                const Reaction<ThermoType>& R = reactions_[i];
                scalar omegaf, omegar;
                R.omega(pi, Ti, c_, celli, omegaf, omegar);

                scalar wf = 0;
                forAll(R.rhs(), s)
                {
                    wf += R.rhs()[s].stoichCoeff*omegaf;
                }
                sumW += wf;
                sumWRateByCTot += sqr(wf);

                scalar wr = 0;
                forAll(R.lhs(), s)
                {
                    wr += R.lhs()[s].stoichCoeff*omegar;
                }
                sumW += wr;
                sumWRateByCTot += sqr(wr);
            }

            tc[celli] =
                sumWRateByCTot == 0 ? vGreat : sumW/sumWRateByCTot*sum(c_);
        }
    }

    ttc.ref().correctBoundaryConditions();

    return ttc;
}


template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::chemistryModel<ThermoType>::Qdot() const
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
Foam::chemistryModel<ThermoType>::calculateRR
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
void Foam::chemistryModel<ThermoType>::calculate()
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
Foam::scalar Foam::chemistryModel<ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    tabulation_.reset();

    const basicSpecieMixture& composition = this->thermo().composition();

    // CPU time analysis
    const clockTime clockTime_ = clockTime();
    clockTime_.timeIncrement();
    scalar solveChemistryCpuTime_ = 0;

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

    // Composition vector (Yi, T, p, deltaT)
    scalarField phiq(nEqns() + 1);
    scalarField Rphiq(nEqns() + 1);

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        scalar pi = p[celli];
        scalar Ti = T[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            c_[i] =
                rhoi*Y_[i].oldTime()[celli]/specieThermos_[i].W();
            c0[i] = c_[i];
            phiq[i] = Y_[i].oldTime()[celli];
        }
        phiq[nSpecie()] = Ti;
        phiq[nSpecie() + 1] = pi;
        phiq[nSpecie() + 2] = deltaT[celli];

        // Initialise time progress
        scalar timeLeft = deltaT[celli];

        // Not sure if this is necessary
        Rphiq = Zero;

        // When tabulation is active (short-circuit evaluation for retrieve)
        // It first tries to retrieve the solution of the system with the
        // information stored through the tabulation method
        if (tabulation_.retrieve(phiq, Rphiq))
        {
            // Retrieved solution stored in Rphiq
            for (label i=0; i<nSpecie(); i++)
            {
                c_[i] = rhoi*Rphiq[i]/specieThermos_[i].W();
            }
        }
        // This position is reached when tabulation is not used OR
        // if the solution is not retrieved.
        // In the latter case, it adds the information to the tabulation
        // (it will either expand the current data or add a new stored point).
        else
        {
            if (mechRedActive_)
            {
                // Reduce mechanism change the number of species (only active)
                mechRed_.reduceMechanism
                (
                    pi,
                    Ti,
                    c_,
                    sc_,
                    cTos_,
                    sToc_,
                    celli
                );
            }

            if (log_)
            {
                // Reset the time
                clockTime_.timeIncrement();
            }

            // Calculate the chemical source terms
            while (timeLeft > small)
            {
                scalar dt = timeLeft;
                if (mechRedActive_)
                {
                    // Solve the reduced set of ODE
                    solve
                    (
                        pi,
                        Ti,
                        sc_,
                        celli,
                        dt,
                        deltaTChem_[celli]
                    );

                    for (label i=0; i<mechRed_.nActiveSpecies(); i++)
                    {
                        c_[sToc_[i]] = sc_[i];
                    }
                }
                else
                {
                    solve
                    (pi, Ti, c_, celli, dt, deltaTChem_[celli]);
                }
                timeLeft -= dt;
            }

            if (log_)
            {
                solveChemistryCpuTime_ += clockTime_.timeIncrement();
            }

            // If tabulation is used, we add the information computed here to
            // the stored points (either expand or add)
            if (tabulation_.tabulates())
            {
                forAll(c_, i)
                {
                    Rphiq[i] = c_[i]/rhoi*specieThermos_[i].W();
                }
                Rphiq[Rphiq.size()-3] = Ti;
                Rphiq[Rphiq.size()-2] = pi;
                Rphiq[Rphiq.size()-1] = deltaT[celli];

                tabulation_.add(phiq, Rphiq, celli, rhoi, deltaT[celli]);
            }

            // When operations are done and if mechanism reduction is active,
            // the number of species (which also affects nEqns) is set back
            // to the total number of species (stored in the mechRed object)
            if (mechRedActive_)
            {
                nSpecie_ = mechRed_.nSpecie();
            }
            deltaTMin = min(deltaTChem_[celli], deltaTMin);

            deltaTChem_[celli] =
                min(deltaTChem_[celli], deltaTChemMax_);
        }

        // Set the RR vector (used in the solver)
        for (label i=0; i<nSpecie_; i++)
        {
            RR_[i][celli] =
                (c_[i] - c0[i])*specieThermos_[i].W()/deltaT[celli];
        }
    }

    if (log_)
    {
        cpuSolveFile_()
            << this->time().userTimeValue()
            << "    " << solveChemistryCpuTime_ << endl;
    }

    mechRed_.update();
    tabulation_.update();

    if (mechRedActive_ && Pstream::parRun())
    {
        List<bool> active(composition.active());
        Pstream::listCombineGather(active, orEqOp<bool>());
        Pstream::listCombineScatter(active);

        forAll(active, i)
        {
            if (active[i])
            {
                composition.setActive(i);
            }
        }
    }

    return deltaTMin;
}


template<class ThermoType>
Foam::scalar Foam::chemistryModel<ThermoType>::solve
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
Foam::scalar Foam::chemistryModel<ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


// ************************************************************************* //
