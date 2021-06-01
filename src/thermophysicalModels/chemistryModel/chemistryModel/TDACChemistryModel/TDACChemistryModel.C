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

#include "TDACChemistryModel.H"
#include "UniformField.H"
#include "localEulerDdtScheme.H"
#include "clockTime.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::TDACChemistryModel<ThermoType>::TDACChemistryModel
(
    const fluidReactionThermo& thermo
)
:
    standardChemistryModel<ThermoType>(thermo),
    timeSteps_(0),
    NsDAC_(this->nSpecie_),
    completeC_(this->nSpecie_, 0),
    reactionsDisabled_(this->reactions_.size(), false),
    completeToSimplifiedIndex_(this->nSpecie_, -1),
    simplifiedToCompleteIndex_(this->nSpecie_),
    tabulationResults_
    (
        IOobject
        (
            thermo.phasePropertyName("TabulationResults"),
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        scalar(0)
    )
{
    const basicSpecieMixture& composition = this->thermo().composition();

    mechRed_ = chemistryReductionMethod<ThermoType>::New
    (
        *this,
        *this
    );

    // When the mechanism reduction method is used, the 'active' flag for every
    // species should be initialised (by default 'active' is true)
    if (mechRed_->active())
    {
        forAll(this->Y(), i)
        {
            IOobject header
            (
                this->Y()[i].name(),
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ
            );

            // Check if the species file is provided, if not set inactive
            // and NO_WRITE
            if (!header.typeHeaderOk<volScalarField>(true))
            {
                composition.setInactive(i);
            }
        }
    }

    tabulation_ = chemistryTabulationMethod<ThermoType>::New
    (
        *this,
        *this
    );

    if (mechRed_->log())
    {
        cpuReduceFile_ = logFile("cpu_reduce.out");
        nActiveSpeciesFile_ = logFile("nActiveSpecies.out");
    }

    if (tabulation_->log())
    {
        cpuAddFile_ = logFile("cpu_add.out");
        cpuGrowFile_ = logFile("cpu_grow.out");
        cpuRetrieveFile_ = logFile("cpu_retrieve.out");
    }

    if (mechRed_->log() || tabulation_->log())
    {
        cpuSolveFile_ = logFile("cpu_solve.out");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::TDACChemistryModel<ThermoType>::~TDACChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::TDACChemistryModel<ThermoType>::omega
(
    const scalar p,
    const scalar T,
    const scalarField& c, // Contains all species even when mechRed is active
    const label li,
    scalarField& dcdt
) const
{
    const bool reduced = mechRed_->active();

    scalar omegaf, omegar;

    dcdt = Zero;

    forAll(this->reactions_, i)
    {
        if (!reactionsDisabled_[i])
        {
            const Reaction<ThermoType>& R = this->reactions_[i];
            const scalar omegaI = R.omega(p, T, c, li, omegaf, omegar);

            forAll(R.lhs(), s)
            {
                const label si =
                    reduced
                  ? completeToSimplifiedIndex_[R.lhs()[s].index]
                  : R.lhs()[s].index;
                const scalar sl = R.lhs()[s].stoichCoeff;
                dcdt[si] -= sl*omegaI;
            }

            forAll(R.rhs(), s)
            {
                const label si =
                    reduced
                  ? completeToSimplifiedIndex_[R.rhs()[s].index]
                  : R.rhs()[s].index;
                const scalar sr = R.rhs()[s].stoichCoeff;
                dcdt[si] += sr*omegaI;
            }
        }
    }
}


template<class ThermoType>
void Foam::TDACChemistryModel<ThermoType>::derivatives
(
    const scalar time,
    const scalarField& c,
    const label li,
    scalarField& dcdt
) const
{
    const bool reduced = mechRed_->active();

    if (reduced)
    {
        // When using DAC, the ODE solver submit a reduced set of species the
        // complete set is used and only the species in the simplified mechanism
        // are updated
        this->c_ = completeC_;

        // Update the concentration of the species in the simplified mechanism
        // the other species remain the same and are used only for third-body
        // efficiencies
        for (label i=0; i<NsDAC_; i++)
        {
            this->c_[simplifiedToCompleteIndex_[i]] = max(c[i], 0);
        }
    }
    else
    {
        for (label i=0; i<this->nSpecie(); i++)
        {
            this->c_[i] = max(c[i], 0);
        }
    }

    const scalar T = c[this->nSpecie_];
    const scalar p = c[this->nSpecie_ + 1];

    dcdt = Zero;

    // Evaluate contributions from reactions
    omega(p, T, this->c_, li, dcdt);

    // Evaluate the effect on the thermodynamic system ...

    // c*Cp
    scalar ccp = 0;
    for (label i=0; i<this->c_.size(); i++)
    {
        ccp += this->c_[i]*this->specieThermos_[i].cp(p, T);
    }

    // dT/dt
    scalar& dTdt = dcdt[this->nSpecie_];
    for (label i=0; i<this->nSpecie_; i++)
    {
        const label si = reduced ? simplifiedToCompleteIndex_[i] : i;
        dTdt -= dcdt[i]*this->specieThermos_[si].ha(p, T);
    }
    dTdt /= ccp;

    // dp/dt = 0 (pressure is assumed constant)
}


template<class ThermoType>
void Foam::TDACChemistryModel<ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    const label li,
    scalarField& dcdt,
    scalarSquareMatrix& J
) const
{
    const bool reduced = mechRed_->active();

    // If the mechanism reduction is active, the computed Jacobian
    // is compact (size of the reduced set of species)
    // but according to the information of the complete set
    // (i.e. for the third-body efficiencies)

    if (reduced)
    {
        this->c_ = completeC_;

        for (label i=0; i<NsDAC_; i++)
        {
            this->c_[simplifiedToCompleteIndex_[i]] = max(c[i], 0);
        }
    }
    else
    {
        forAll(this->c_, i)
        {
            this->c_[i] = max(c[i], 0);
        }
    }

    const scalar T = c[this->nSpecie_];
    const scalar p = c[this->nSpecie_ + 1];

    dcdt = Zero;
    J = Zero;

    // Evaluate contributions from reactions
    forAll(this->reactions_, ri)
    {
        if (!reactionsDisabled_[ri])
        {
            const Reaction<ThermoType>& R = this->reactions_[ri];
            scalar omegaI, kfwd, kbwd;
            R.dwdc
            (
                p,
                T,
                this->c_,
                li,
                J,
                dcdt,
                omegaI,
                kfwd,
                kbwd,
                reduced,
                completeToSimplifiedIndex_
            );
            R.dwdT
            (
                p,
                T,
                this->c_,
                li,
                omegaI,
                kfwd,
                kbwd,
                J,
                reduced,
                completeToSimplifiedIndex_,
                this->nSpecie_
            );
        }
    }

    // Evaluate the effect on the thermodynamic system ...

    // c*Cp
    scalar ccp = 0, dccpdT = 0;
    forAll(this->c_, i)
    {
        ccp += this->c_[i]*this->specieThermos_[i].cp(p, T);
        dccpdT += this->c_[i]*this->specieThermos_[i].dcpdT(p, T);
    }

    // dT/dt
    scalar& dTdt = dcdt[this->nSpecie_];
    for (label i=0; i<this->nSpecie_; i++)
    {
        const label si = reduced ? simplifiedToCompleteIndex_[i] : i;
        dTdt -= dcdt[i]*this->specieThermos_[si].ha(p, T);
    }
    dTdt /= ccp;

    // dp/dt = 0 (pressure is assumed constant)

    // d(dTdt)/dc
    for (label i = 0; i < this->nSpecie_; i++)
    {
        scalar& d2Tdtdci = J(this->nSpecie_, i);
        for (label j = 0; j < this->nSpecie_; j++)
        {
            const scalar d2cjdtdci = J(j, i);
            const label sj = reduced ? simplifiedToCompleteIndex_[j] : j;
            d2Tdtdci -= d2cjdtdci*this->specieThermos_[sj].ha(p, T);
        }
        const label si = reduced ? simplifiedToCompleteIndex_[i] : i;
        d2Tdtdci -= this->specieThermos_[si].cp(p, T)*dTdt;
        d2Tdtdci /= ccp;
    }

    // d(dTdt)/dT
    scalar& d2TdtdT = J(this->nSpecie_, this->nSpecie_);
    for (label i = 0; i < this->nSpecie_; i++)
    {
        const scalar d2cidtdT = J(i, this->nSpecie_);
        const label si = reduced ? simplifiedToCompleteIndex_[i] : i;
        d2TdtdT -=
            dcdt[i]*this->specieThermos_[si].cp(p, T)
          + d2cidtdT*this->specieThermos_[si].ha(p, T);
    }
    d2TdtdT -= dTdt*dccpdT;
    d2TdtdT /= ccp;

    // d(dpdt)/dc = 0 (pressure is assumed constant)

    // d(dpdt)/dT = 0 (pressure is assumed constant)
}


template<class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::TDACChemistryModel<ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    // Increment counter of time-step
    timeSteps_++;

    const bool reduced = mechRed_->active();

    const basicSpecieMixture& composition = this->thermo().composition();

    // CPU time analysis
    const clockTime clockTime_= clockTime();
    clockTime_.timeIncrement();
    scalar reduceMechCpuTime_ = 0;
    scalar addNewLeafCpuTime_ = 0;
    scalar growCpuTime_ = 0;
    scalar solveChemistryCpuTime_ = 0;
    scalar searchISATCpuTime_ = 0;

    this->resetTabulationResults();

    // Average number of active species
    scalar nActiveSpecies = 0;
    scalar nAvg = 0;

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

    scalarField c(this->nSpecie_);
    scalarField c0(this->nSpecie_);

    // Composition vector (Yi, T, p, deltaT)
    scalarField phiq(this->nEqns() + 1);
    scalarField Rphiq(this->nEqns() + 1);

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        scalar pi = p[celli];
        scalar Ti = T[celli];

        for (label i=0; i<this->nSpecie_; i++)
        {
            c[i] = c0[i] =
                rhoi*this->Y_[i].oldTime()[celli]/this->specieThermos_[i].W();
            phiq[i] = this->Y()[i].oldTime()[celli];
        }
        phiq[this->nSpecie()] = Ti;
        phiq[this->nSpecie() + 1] = pi;
        phiq[this->nSpecie() + 2] = deltaT[celli];

        // Initialise time progress
        scalar timeLeft = deltaT[celli];

        // Not sure if this is necessary
        Rphiq = Zero;

        clockTime_.timeIncrement();

        // When tabulation is active (short-circuit evaluation for retrieve)
        // It first tries to retrieve the solution of the system with the
        // information stored through the tabulation method
        if (tabulation_->active() && tabulation_->retrieve(phiq, Rphiq))
        {
            // Retrieved solution stored in Rphiq
            for (label i=0; i<this->nSpecie(); i++)
            {
                c[i] = rhoi*Rphiq[i]/this->specieThermos_[i].W();
            }

            searchISATCpuTime_ += clockTime_.timeIncrement();
        }
        // This position is reached when tabulation is not used OR
        // if the solution is not retrieved.
        // In the latter case, it adds the information to the tabulation
        // (it will either expand the current data or add a new stored point).
        else
        {
            // Reset the time
            clockTime_.timeIncrement();

            if (reduced)
            {
                // Reduce mechanism change the number of species (only active)
                mechRed_->reduceMechanism(pi, Ti, c, celli);
                nActiveSpecies += mechRed_->NsSimp();
                nAvg++;
                reduceMechCpuTime_ += clockTime_.timeIncrement();
            }

            // Calculate the chemical source terms
            while (timeLeft > small)
            {
                scalar dt = timeLeft;
                if (reduced)
                {
                    // completeC_ used in the overridden ODE methods
                    // to update only the active species
                    completeC_ = c;

                    // Solve the reduced set of ODE
                    this->solve
                    (
                        pi,
                        Ti,
                        simplifiedC_,
                        celli,
                        dt,
                        this->deltaTChem_[celli]
                    );

                    for (label i=0; i<NsDAC_; i++)
                    {
                        c[simplifiedToCompleteIndex_[i]] = simplifiedC_[i];
                    }
                }
                else
                {
                    this->solve(pi, Ti, c, celli, dt, this->deltaTChem_[celli]);
                }
                timeLeft -= dt;
            }

            {
                solveChemistryCpuTime_ += clockTime_.timeIncrement();
            }

            // If tabulation is used, we add the information computed here to
            // the stored points (either expand or add)
            if (tabulation_->active())
            {
                forAll(c, i)
                {
                    Rphiq[i] = c[i]/rhoi*this->specieThermos_[i].W();
                }
                Rphiq[Rphiq.size()-3] = Ti;
                Rphiq[Rphiq.size()-2] = pi;
                Rphiq[Rphiq.size()-1] = deltaT[celli];

                label growOrAdd =
                    tabulation_->add(phiq, Rphiq, celli, rhoi, deltaT[celli]);
                if (growOrAdd)
                {
                    this->setTabulationResultsAdd(celli);
                    addNewLeafCpuTime_ += clockTime_.timeIncrement();
                }
                else
                {
                    this->setTabulationResultsGrow(celli);
                    growCpuTime_ += clockTime_.timeIncrement();
                }
            }

            // When operations are done and if mechanism reduction is active,
            // the number of species (which also affects nEqns) is set back
            // to the total number of species (stored in the mechRed object)
            if (reduced)
            {
                this->nSpecie_ = mechRed_->nSpecie();
            }
            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

            this->deltaTChem_[celli] =
                min(this->deltaTChem_[celli], this->deltaTChemMax_);
        }

        // Set the RR vector (used in the solver)
        for (label i=0; i<this->nSpecie_; i++)
        {
            this->RR_[i][celli] =
                (c[i] - c0[i])*this->specieThermos_[i].W()/deltaT[celli];
        }
    }

    if (mechRed_->log() || tabulation_->log())
    {
        cpuSolveFile_()
            << this->time().timeOutputValue()
            << "    " << solveChemistryCpuTime_ << endl;
    }

    if (mechRed_->log())
    {
        cpuReduceFile_()
            << this->time().timeOutputValue()
            << "    " << reduceMechCpuTime_ << endl;
    }

    if (tabulation_->active())
    {
        // Every time-step, look if the tabulation should be updated
        tabulation_->update();

        // Write the performance of the tabulation
        tabulation_->writePerformance();

        if (tabulation_->log())
        {
            cpuRetrieveFile_()
                << this->time().timeOutputValue()
                << "    " << searchISATCpuTime_ << endl;

            cpuGrowFile_()
                << this->time().timeOutputValue()
                << "    " << growCpuTime_ << endl;

            cpuAddFile_()
                << this->time().timeOutputValue()
                << "    " << addNewLeafCpuTime_ << endl;
        }
    }

    if (reduced && nAvg && mechRed_->log())
    {
        // Write average number of species
        nActiveSpeciesFile_()
            << this->time().timeOutputValue()
            << "    " << nActiveSpecies/nAvg << endl;
    }

    if (Pstream::parRun())
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
Foam::scalar Foam::TDACChemistryModel<ThermoType>::solve
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
Foam::scalar Foam::TDACChemistryModel<ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template<class ThermoType>
void Foam::TDACChemistryModel<ThermoType>::
setTabulationResultsAdd
(
    const label celli
)
{
    tabulationResults_[celli] = 0;
}


template<class ThermoType>
void Foam::TDACChemistryModel<ThermoType>::
setTabulationResultsGrow(const label celli)
{
    tabulationResults_[celli] = 1;
}


template<class ThermoType>
void Foam::TDACChemistryModel<ThermoType>::
setTabulationResultsRetrieve
(
    const label celli
)
{
    tabulationResults_[celli] = 2;
}


// ************************************************************************* //
