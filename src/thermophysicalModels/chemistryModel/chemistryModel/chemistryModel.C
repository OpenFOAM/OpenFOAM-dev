/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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
#include "cpuLoad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::chemistryModel<ThermoType>::chemistryModel
(
    const fluidMulticomponentThermo& thermo
)
:
    odeChemistryModel(thermo),
    log_(this->lookupOrDefault("log", false)),
    loadBalancing_(this->lookupOrDefault("loadBalancing", false)),
    jacobianType_
    (
        this->found("jacobian")
      ? jacobianTypeNames_.read(this->lookup("jacobian"))
      : jacobianType::fast
    ),
    mixture_(refCast<const multicomponentMixture<ThermoType>>(this->thermo())),
    specieThermos_(mixture_.specieThermos()),
    reactions_(mixture_.species(), specieThermos_, this->mesh(), *this),
    RR_(nSpecie_),
    Y_(nSpecie_),
    c_(nSpecie_),
    YTpWork_(scalarField(nSpecie_ + 2)),
    YTpYTpWork_(scalarSquareMatrix(nSpecie_ + 2)),
    mechRedPtr_
    (
        chemistryReductionMethod<ThermoType>::New
        (
            *this,
            *this
        )
    ),
    mechRed_(*mechRedPtr_),
    tabulationPtr_(chemistryTabulationMethod::New(*this, *this)),
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
                    "RR." + Yvf_[fieldi].name(),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermo.mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0)
            )
        );
    }

    Info<< "chemistryModel: Number of species = " << nSpecie_
        << " and reactions = " << nReaction() << endl;

    // When the mechanism reduction method is used, the 'active' flag for every
    // species should be initialised (by default 'active' is true)
    if (reduction_)
    {
        const basicSpecieMixture& composition = this->thermo().composition();

        forAll(Yvf_, i)
        {
            typeIOobject<volScalarField> header
            (
                Yvf_[i].name(),
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
void Foam::chemistryModel<ThermoType>::derivatives
(
    const scalar time,
    const scalarField& YTp,
    const label li,
    scalarField& dYTpdt
) const
{
    if (reduction_)
    {
        forAll(sToc_, i)
        {
            Y_[sToc_[i]] = max(YTp[i], 0);
        }
    }
    else
    {
        forAll(Y_, i)
        {
            Y_[i] = max(YTp[i], 0);
        }
    }

    const scalar T = YTp[nSpecie_];
    const scalar p = YTp[nSpecie_ + 1];

    // Evaluate the mixture density
    scalar rhoM = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        rhoM += Y_[i]/specieThermos_[i].rho(p, T);
    }
    rhoM = 1/rhoM;

    // Evaluate the concentrations
    for (label i=0; i<Y_.size(); i ++)
    {
        c_[i] = rhoM/specieThermos_[i].W()*Y_[i];
    }

    // Evaluate contributions from reactions
    dYTpdt = Zero;
    forAll(reactions_, ri)
    {
        if (!mechRed_.reactionDisabled(ri))
        {
            reactions_[ri].dNdtByV
            (
                p,
                T,
                c_,
                li,
                dYTpdt,
                reduction_,
                cTos_,
                0
            );
        }
    }

    // Reactions return dNdtByV, so we need to convert the result to dYdt
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar WiByrhoM = specieThermos_[sToc(i)].W()/rhoM;
        scalar& dYidt = dYTpdt[i];
        dYidt *= WiByrhoM;
    }

    // Evaluate the effect on the thermodynamic system ...

    // Evaluate the mixture Cp
    scalar CpM = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        CpM += Y_[i]*specieThermos_[i].Cp(p, T);
    }

    // dT/dt
    scalar& dTdt = dYTpdt[nSpecie_];
    for (label i=0; i<nSpecie_; i++)
    {
        dTdt -= dYTpdt[i]*specieThermos_[sToc(i)].Ha(p, T);
    }
    dTdt /= CpM;

    // dp/dt = 0 (pressure is assumed constant)
    scalar& dpdt = dYTpdt[nSpecie_ + 1];
    dpdt = 0;
}


template<class ThermoType>
void Foam::chemistryModel<ThermoType>::jacobian
(
    const scalar t,
    const scalarField& YTp,
    const label li,
    scalarField& dYTpdt,
    scalarSquareMatrix& J
) const
{
    if (reduction_)
    {
        forAll(sToc_, i)
        {
            Y_[sToc_[i]] = max(YTp[i], 0);
        }
    }
    else
    {
        forAll(c_, i)
        {
            Y_[i] = max(YTp[i], 0);
        }
    }

    const scalar T = YTp[nSpecie_];
    const scalar p = YTp[nSpecie_ + 1];

    // Evaluate the specific volumes and mixture density
    scalarField& v = YTpWork_[0];
    for (label i=0; i<Y_.size(); i++)
    {
        v[i] = 1/specieThermos_[i].rho(p, T);
    }
    scalar rhoM = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        rhoM += Y_[i]*v[i];
    }
    rhoM = 1/rhoM;

    // Evaluate the concentrations
    for (label i=0; i<Y_.size(); i ++)
    {
        c_[i] = rhoM/specieThermos_[i].W()*Y_[i];
    }

    // Evaluate the derivatives of concentration w.r.t. mass fraction
    scalarSquareMatrix& dcdY = YTpYTpWork_[0];
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar rhoMByWi = rhoM/specieThermos_[sToc(i)].W();
        switch (jacobianType_)
        {
            case jacobianType::fast:
                {
                    dcdY(i, i) = rhoMByWi;
                }
                break;
            case jacobianType::exact:
                for (label j=0; j<nSpecie_; j++)
                {
                    dcdY(i, j) =
                        rhoMByWi*((i == j) - rhoM*v[sToc(j)]*Y_[sToc(i)]);
                }
                break;
        }
    }

    // Evaluate the mixture thermal expansion coefficient
    scalar alphavM = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        alphavM += Y_[i]*rhoM*v[i]*specieThermos_[i].alphav(p, T);
    }

    // Evaluate contributions from reactions
    dYTpdt = Zero;
    scalarSquareMatrix& ddNdtByVdcTp = YTpYTpWork_[1];
    for (label i=0; i<nSpecie_ + 2; i++)
    {
        for (label j=0; j<nSpecie_ + 2; j++)
        {
            ddNdtByVdcTp[i][j] = 0;
        }
    }
    forAll(reactions_, ri)
    {
        if (!mechRed_.reactionDisabled(ri))
        {
            reactions_[ri].ddNdtByVdcTp
            (
                p,
                T,
                c_,
                li,
                dYTpdt,
                ddNdtByVdcTp,
                reduction_,
                cTos_,
                0,
                nSpecie_,
                YTpWork_[1],
                YTpWork_[2]
            );
        }
    }

    // Reactions return dNdtByV, so we need to convert the result to dYdt
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar WiByrhoM = specieThermos_[sToc(i)].W()/rhoM;
        scalar& dYidt = dYTpdt[i];
        dYidt *= WiByrhoM;

        for (label j=0; j<nSpecie_; j++)
        {
            scalar ddNidtByVdYj = 0;
            switch (jacobianType_)
            {
                case jacobianType::fast:
                    {
                        const scalar ddNidtByVdcj = ddNdtByVdcTp(i, j);
                        ddNidtByVdYj = ddNidtByVdcj*dcdY(j, j);
                    }
                    break;
                case jacobianType::exact:
                    for (label k=0; k<nSpecie_; k++)
                    {
                        const scalar ddNidtByVdck = ddNdtByVdcTp(i, k);
                        ddNidtByVdYj += ddNidtByVdck*dcdY(k, j);
                    }
                    break;
            }

            scalar& ddYidtdYj = J(i, j);
            ddYidtdYj = WiByrhoM*ddNidtByVdYj + rhoM*v[sToc(j)]*dYidt;
        }

        scalar ddNidtByVdT = ddNdtByVdcTp(i, nSpecie_);
        for (label j=0; j<nSpecie_; j++)
        {
            const scalar ddNidtByVdcj = ddNdtByVdcTp(i, j);
            ddNidtByVdT -= ddNidtByVdcj*c_[sToc(j)]*alphavM;
        }

        scalar& ddYidtdT = J(i, nSpecie_);
        ddYidtdT = WiByrhoM*ddNidtByVdT + alphavM*dYidt;

        scalar& ddYidtdp = J(i, nSpecie_ + 1);
        ddYidtdp = 0;
    }

    // Evaluate the effect on the thermodynamic system ...

    // Evaluate the mixture Cp and its derivative
    scalarField& Cp = YTpWork_[3];
    scalar CpM = 0, dCpMdT = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        Cp[i] = specieThermos_[i].Cp(p, T);
        CpM += Y_[i]*Cp[i];
        dCpMdT += Y_[i]*specieThermos_[i].dCpdT(p, T);
    }

    // dT/dt
    scalarField& Ha = YTpWork_[4];
    scalar& dTdt = dYTpdt[nSpecie_];
    for (label i=0; i<nSpecie_; i++)
    {
        Ha[sToc(i)] = specieThermos_[sToc(i)].Ha(p, T);
        dTdt -= dYTpdt[i]*Ha[sToc(i)];
    }
    dTdt /= CpM;

    // dp/dt = 0 (pressure is assumed constant)
    scalar& dpdt = dYTpdt[nSpecie_ + 1];
    dpdt = 0;

    // d(dTdt)/dY
    for (label i=0; i<nSpecie_; i++)
    {
        scalar& ddTdtdYi = J(nSpecie_, i);
        ddTdtdYi = 0;
        for (label j=0; j<nSpecie_; j++)
        {
            const scalar ddYjdtdYi = J(j, i);
            ddTdtdYi -= ddYjdtdYi*Ha[sToc(j)];
        }
        ddTdtdYi -= Cp[sToc(i)]*dTdt;
        ddTdtdYi /= CpM;
    }

    // d(dTdt)/dT
    scalar& ddTdtdT = J(nSpecie_, nSpecie_);
    ddTdtdT = 0;
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar dYidt = dYTpdt[i];
        const scalar ddYidtdT = J(i, nSpecie_);
        ddTdtdT -= dYidt*Cp[sToc(i)] + ddYidtdT*Ha[sToc(i)];
    }
    ddTdtdT -= dTdt*dCpMdT;
    ddTdtdT /= CpM;

    // d(dTdt)/dp = 0 (pressure is assumed constant)
    scalar& ddTdtdp = J(nSpecie_, nSpecie_ + 1);
    ddTdtdp = 0;

    // d(dpdt)/dYiTp = 0 (pressure is assumed constant)
    for (label i=0; i<nSpecie_ + 2; i++)
    {
        scalar& ddpdtdYiTp = J(nSpecie_ + 1, i);
        ddpdtdYiTp = 0;
    }
}


template<class ThermoType>
Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::chemistryModel<ThermoType>::reactionRR
(
    const label reactioni
) const
{
    PtrList<volScalarField::Internal> RR(nSpecie_);
    for (label i=0; i<nSpecie_; i++)
    {
        RR.set
        (
            i,
            volScalarField::Internal::New
            (
                "RR." + Yvf_[i].name(),
                this->mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0)
            ).ptr()
        );
    }

    if (!this->chemistry_ || mechRed_.reactionDisabled(reactioni))
    {
        return RR;
    }

    tmp<volScalarField> trhovf(this->thermo().rho());
    const volScalarField& rhovf = trhovf();

    const volScalarField& Tvf = this->thermo().T();
    const volScalarField& pvf = this->thermo().p();

    scalarField& dNdtByV = YTpWork_[0];

    reactionEvaluationScope scope(*this);

    const Reaction<ThermoType>& R = reactions_[reactioni];

    forAll(rhovf, celli)
    {
        const scalar rho = rhovf[celli];
        const scalar T = Tvf[celli];
        const scalar p = pvf[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Yvf_[i][celli];
            c_[i] = rho*Yi/specieThermos_[i].W();
        }

        dNdtByV = Zero;

        R.dNdtByV
        (
            p,
            T,
            c_,
            celli,
            dNdtByV,
            reduction_,
            cTos_,
            0
        );

        for (label i=0; i<nSpecie_; i++)
        {
            RR[i][celli] = dNdtByV[i]*specieThermos_[i].W();
        }
    }

    return RR;
}


template<class ThermoType>
void Foam::chemistryModel<ThermoType>::calculate()
{
    if (!this->chemistry_)
    {
        return;
    }

    tmp<volScalarField> trhovf(this->thermo().rho());
    const volScalarField& rhovf = trhovf();

    const volScalarField& Tvf = this->thermo().T();
    const volScalarField& pvf = this->thermo().p();

    scalarField& dNdtByV = YTpWork_[0];

    reactionEvaluationScope scope(*this);

    forAll(rhovf, celli)
    {
        const scalar rho = rhovf[celli];
        const scalar T = Tvf[celli];
        const scalar p = pvf[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Yvf_[i][celli];
            c_[i] = rho*Yi/specieThermos_[i].W();
        }

        dNdtByV = Zero;

        forAll(reactions_, ri)
        {
            if (!mechRed_.reactionDisabled(ri))
            {
                reactions_[ri].dNdtByV
                (
                    p,
                    T,
                    c_,
                    celli,
                    dNdtByV,
                    reduction_,
                    cTos_,
                    0
                );
            }
        }

        for (label i=0; i<mechRed_.nActiveSpecies(); i++)
        {
            RR_[sToc(i)][celli] = dNdtByV[i]*specieThermos_[sToc(i)].W();
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
    optionalCpuLoad& chemistryCpuTime
    (
        optionalCpuLoad::New(this->mesh(), "chemistryCpuTime", loadBalancing_)
    );

    // CPU time logging
    cpuTime solveCpuTime;
    scalar totalSolveCpuTime = 0;

    if (!this->chemistry_)
    {
        return great;
    }

    const volScalarField& rho0vf =
        this->mesh().template lookupObject<volScalarField>
        (
            this->thermo().phasePropertyName("rho")
        ).oldTime();

    const volScalarField& T0vf = this->thermo().T().oldTime();
    const volScalarField& p0vf = this->thermo().p().oldTime();

    reactionEvaluationScope scope(*this);

    scalarField Y0(nSpecie_);

    // Composition vector (Yi, T, p, deltaT)
    scalarField phiq(nEqns() + 1);
    scalarField Rphiq(nEqns() + 1);

    // Minimum chemical timestep
    scalar deltaTMin = great;

    tabulation_.reset();
    chemistryCpuTime.reset();

    forAll(rho0vf, celli)
    {
        const scalar rho0 = rho0vf[celli];

        scalar p = p0vf[celli];
        scalar T = T0vf[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            Y_[i] = Y0[i] = Yvf_[i].oldTime()[celli];
        }

        for (label i=0; i<nSpecie_; i++)
        {
            phiq[i] = Yvf_[i].oldTime()[celli];
        }
        phiq[nSpecie()] = T;
        phiq[nSpecie() + 1] = p;
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
                Y_[i] = Rphiq[i];
            }
            T = Rphiq[nSpecie()];
            p = Rphiq[nSpecie() + 1];
        }
        // This position is reached when tabulation is not used OR
        // if the solution is not retrieved.
        // In the latter case, it adds the information to the tabulation
        // (it will either expand the current data or add a new stored point).
        else
        {
            if (reduction_)
            {
                // Compute concentrations
                for (label i=0; i<nSpecie_; i++)
                {
                    c_[i] = rho0*Y_[i]/specieThermos_[i].W();
                }

                // Reduce mechanism change the number of species (only active)
                mechRed_.reduceMechanism(p, T, c_, cTos_, sToc_, celli);

                // Set the simplified mass fraction field
                sY_.setSize(nSpecie_);
                for (label i=0; i<nSpecie_; i++)
                {
                    sY_[i] = Y_[sToc(i)];
                }
            }

            if (log_)
            {
                // Reset the solve time
                solveCpuTime.cpuTimeIncrement();
            }

            // Calculate the chemical source terms
            while (timeLeft > small)
            {
                scalar dt = timeLeft;
                if (reduction_)
                {
                    // Solve the reduced set of ODE
                    solve
                    (
                        p,
                        T,
                        sY_,
                        celli,
                        dt,
                        deltaTChem_[celli]
                    );

                    for (label i=0; i<mechRed_.nActiveSpecies(); i++)
                    {
                        Y_[sToc_[i]] = sY_[i];
                    }
                }
                else
                {
                    solve(p, T, Y_, celli, dt, deltaTChem_[celli]);
                }
                timeLeft -= dt;
            }

            if (log_)
            {
                totalSolveCpuTime += solveCpuTime.cpuTimeIncrement();
            }

            // If tabulation is used, we add the information computed here to
            // the stored points (either expand or add)
            if (tabulation_.tabulates())
            {
                forAll(Y_, i)
                {
                    Rphiq[i] = Y_[i];
                }
                Rphiq[Rphiq.size()-3] = T;
                Rphiq[Rphiq.size()-2] = p;
                Rphiq[Rphiq.size()-1] = deltaT[celli];

                tabulation_.add
                (
                    phiq,
                    Rphiq,
                    mechRed_.nActiveSpecies(),
                    celli,
                    deltaT[celli]
                );
            }

            // When operations are done and if mechanism reduction is active,
            // the number of species (which also affects nEqns) is set back
            // to the total number of species (stored in the mechRed object)
            if (reduction_)
            {
                setNSpecie(mechRed_.nSpecie());
            }

            deltaTMin = min(deltaTChem_[celli], deltaTMin);
            deltaTChem_[celli] = min(deltaTChem_[celli], deltaTChemMax_);
        }

        // Set the RR vector (used in the solver)
        for (label i=0; i<nSpecie_; i++)
        {
            RR_[i][celli] = rho0*(Y_[i] - Y0[i])/deltaT[celli];
        }

        if (loadBalancing_)
        {
            chemistryCpuTime.cpuTimeIncrement(celli);
        }
    }

    if (log_)
    {
        cpuSolveFile_()
            << this->time().userTimeValue()
            << "    " << totalSolveCpuTime << endl;
    }

    mechRed_.update();
    tabulation_.update();

    if (reduction_ && Pstream::parRun())
    {
        const basicSpecieMixture& composition = this->thermo().composition();

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

    if (!this->chemistry_)
    {
        ttc.ref().correctBoundaryConditions();
        return ttc;
    }

    tmp<volScalarField> trhovf(this->thermo().rho());
    const volScalarField& rhovf = trhovf();

    const volScalarField& Tvf = this->thermo().T();
    const volScalarField& pvf = this->thermo().p();

    reactionEvaluationScope scope(*this);

    forAll(rhovf, celli)
    {
        const scalar rho = rhovf[celli];
        const scalar T = Tvf[celli];
        const scalar p = pvf[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            c_[i] = rho*Yvf_[i][celli]/specieThermos_[i].W();
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
            R.omega(p, T, c_, celli, omegaf, omegar);

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

    if (!this->chemistry_)
    {
        return tQdot;
    }

    reactionEvaluationScope scope(*this);

    scalarField& Qdot = tQdot.ref();

    forAll(Yvf_, i)
    {
        forAll(Qdot, celli)
        {
            const scalar hi = specieThermos_[i].Hf();
            Qdot[celli] -= hi*RR_[i][celli];
        }
    }

    return tQdot;
}


// ************************************************************************* //
