/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2026 OpenFOAM Foundation
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

#include "standard_chemistryModel.H"
#include "chemistryReductionMethod.H"
#include "chemistryTabulationMethod.H"
#include "cpuLoad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace chemistryModels
{
    const Foam::NamedEnum<standard::jacobianType, 2>
    standard::jacobianTypeNames
    {
        "fast",
        "exact"
    };

    defineTypeNameAndDebug(standard, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalarList Foam::chemistryModels::standard::specieWs
(
    const fluidMulticomponentThermo& thermo
)
{
    scalarList Ws(thermo.species().size());
    forAll(Ws, speciei)
    {
        Ws[speciei] = thermo.WiValue(speciei);
    }
    return Ws;
}


Foam::scalarList Foam::chemistryModels::standard::specieHfs
(
    const fluidMulticomponentThermo& thermo
)
{
    scalarList hfs(thermo.species().size());
    forAll(hfs, speciei)
    {
        hfs[speciei] = thermo.hfiValue(speciei);
    }
    return hfs;
}


template<class DeltaTType>
Foam::scalar Foam::chemistryModels::standard::solve
(
    const DeltaTType& deltaT
)
{
    optionalCpuLoad& chemistryCpuLoad
    (
        optionalCpuLoad::New(name() + ":cpuLoad", this->mesh(), cpuLoad_)
    );

    // CPU time logging
    cpuTime solveCpuTime;
    scalar totalSolveCpuTime = 0;

    if (!this->chemistry_)
    {
        return great;
    }

    if (!zone_.all())
    {
        forAll(RR_, fieldi)
        {
            RR_[fieldi] = Zero;
        }
    }

    const volScalarField& rho0vf =
        this->mesh().template lookupObject<volScalarField>
        (
            this->thermo().phasePropertyName("rho")
        ).oldTime();

    const volScalarField& T0vf = this->thermo().T().oldTime();
    const volScalarField& p0vf = this->thermo().p().oldTime();

    reactionEvaluationScope scope(*this);

    // Tabulated state vector (Yi, T, p, deltaT)
    scalarField phiq;
    scalarField Rphiq;
    if (tabulation_.tabulates())
    {
        phiq.setSize(nEqns() + 2);
        Rphiq.setSize(nEqns() + 2);
    }

    // Minimum chemical timestep
    scalar deltaTMin = great;

    tabulation_.reset();
    chemistryCpuLoad.resetCpuTime();

    zone_.regenerate();
    const label nZoneCells = zone_.nCells();
    for(label zci = 0; zci<nZoneCells; zci++)
    {
        const label celli = zone_.celli(zci);

        const scalar rho0 = rho0vf[celli];

        p_ = p0vf[celli];
        scalar T = T0vf[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            Y_[i] = Yvf_[i].oldTime()[celli];
        }

        if (tabulation_.tabulates())
        {
            for (label i=0; i<nSpecie_; i++)
            {
                phiq[i] = Yvf_[i].oldTime()[celli];
            }
            phiq[nSpecie()] = T;
            phiq[nSpecie() + 1] = p_;
            phiq[nSpecie() + 2] = deltaT[celli];
        }

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
            // and normalise the mass-fractions
            scalar Yt = 0;
            for (label i=0; i<nSpecie(); i++)
            {
                Yt += Rphiq[i];
            }
            for (label i=0; i<nSpecie(); i++)
            {
                Y_[i] = Rphiq[i]/Yt;
            }
            T = Rphiq[nSpecie()];
        }
        // This position is reached when tabulation is not used OR
        // if the solution is not retrieved.
        // In the latter case, it adds the information to the tabulation
        // (it will either expand the current data or add a new stored point).
        else
        {
            if (reduced_)
            {
                // Compute concentrations
                for (label i=0; i<nSpecie_; i++)
                {
                    s_.c[i] = rho0*Y_[i]/specieWs_[i];
                }

                s_.reset(p_, T, celli);

                // Reduce mechanism change the number of species (only active)
                reduction_.reduceMechanism(s_, c2s_, sToc_);

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
                if (reduced_)
                {
                    // Solve the reduced set of ODE
                    solve
                    (
                        p_,
                        T,
                        sY_,
                        celli,
                        dt,
                        deltaTChem_[celli]
                    );

                    for (label i=0; i<reduction_.nActiveSpecies(); i++)
                    {
                        Y_[sToc_[i]] = sY_[i];
                    }
                }
                else
                {
                    solve(p_, T, Y_, celli, dt, deltaTChem_[celli]);
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
                Rphiq[Rphiq.size()-2] = p_;
                Rphiq[Rphiq.size()-1] = deltaT[celli];

                tabulation_.add
                (
                    phiq,
                    Rphiq,
                    reduction_.nActiveSpecies(),
                    celli,
                    deltaT[celli]
                );
            }

            // When operations are done and if mechanism reduction is active,
            // the number of species (which also affects nEqns) is set back
            // to the total number of species (stored in the reduction object)
            if (reduced_)
            {
                setNSpecie(reduction_.nSpecie());
            }

            deltaTMin = min(deltaTChem_[celli], deltaTMin);
            deltaTChem_[celli] = min(deltaTChem_[celli], deltaTChemMax_);
        }

        // Set the RR vector (used in the solver)
        for (label i=0; i<nSpecie_; i++)
        {
            RR_[i][celli] =
                rho0*(Y_[i] - Yvf_[i].oldTime()[celli])/deltaT[celli];
        }

        if (cpuLoad_)
        {
            chemistryCpuLoad.cpuTimeIncrement(celli);
        }
    }

    if (log_)
    {
        cpuSolveFile_()
            << this->time().userTimeValue()
            << "    " << totalSolveCpuTime << endl;
    }

    reduction_.update();
    tabulation_.update();

    if (reduced_)
    {
        this->thermo().syncActiveSpecies();
    }

    return deltaTMin;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chemistryModels::standard::standard
(
    const fluidMulticomponentThermo& thermo,
    const specieThermosRef& specieThermos
)
:
    chemistryModel(thermo),
    ODESystem(),
    Yvf_(this->thermo().Y()),
    nSpecie_(Yvf_.size()),
    reactions_(thermo.species(), specieThermos, this->mesh(), *this),
    reductionPtr_(chemistryReductionMethod::New(*this, *this)),
    reduction_(*reductionPtr_),
    reduced_(reduction_.active()),
    tabulationPtr_(chemistryTabulationMethod::New(*this, *this)),
    tabulation_(*tabulationPtr_),
    zone_(mesh(), *this, true, true),
    specieWs_(specieWs(thermo)),
    specieHfs_(specieHfs(thermo)),
    c2s_(nSpecie_, -1),
    sToc_(nSpecie_),
    odeSolver_(ODESolver::New(*this, typeDict("ode"))),
    cT_(nEqns()),
    p_(0),
    RR_(nSpecie_),
    Y_(nSpecie_),
    sY_(nSpecie_),
    s_(nSpecie_),
    YTWork_(scalarField(nSpecie_ + 1)),
    YTYTWork_(scalarSquareMatrix(nSpecie_ + 1)),
    log_(this->lookupOrDefault("log", false)),
    cpuLoad_(this->lookupOrDefault("cpuLoad", false))
{
    // Create the fields for the chemistry sources
    forAll(RR_, fieldi)
    {
        RR_.set
        (
            fieldi,
            new volInternalScalarField
            (
                IOobject
                (
                    "RR." + Yvf_[fieldi].name(),
                    this->mesh().time().name(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermo.mesh(),
                dimensionedScalar(dimensions::density/dimensions::time, 0)
            )
        );
    }

    Info<< "chemistryModel: Number of species = " << nSpecie_
        << " and reactions = " << nReaction() << endl;

    // When the mechanism reduction method is used, the 'active' flag for every
    // species should be initialised (by default 'active' is true)
    if (reduced_)
    {
        forAll(Yvf_, i)
        {
            typeIOobject<volScalarField> header
            (
                Yvf_[i].name(),
                this->mesh().time().name(),
                this->mesh(),
                IOobject::NO_READ
            );

            // Check if the species file is provided, if not set inactive
            // and NO_WRITE
            if (!header.headerOk())
            {
                this->thermo().setSpecieInactive(i);
            }
        }

        this->thermo().syncActiveSpecies();
    }

    if (log_)
    {
        cpuSolveFile_ = logFile("cpu_solve.out");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::chemistryModels::standard::~standard()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::fvMesh>>
Foam::chemistryModels::standard::reactionRR
(
    const label reactioni
) const
{
    tmp<volInternalScalarField> tRR =
        volInternalScalarField::New
        (
            "RR:" + reactions_[reactioni].name(),
            this->mesh(),
            dimensionedScalar
            (
                dimensions::moles/dimensions::volume/dimensions::time,
                0
            )
        );
    volInternalScalarField& RR = tRR.ref();

    if (!this->chemistry_ || reduction_.reactionDisabled(reactioni))
    {
        return tRR;
    }

    tmp<volScalarField> trhovf(this->thermo().rho());
    const volScalarField& rhovf = trhovf();

    const volScalarField& Tvf = this->thermo().T();
    const volScalarField& pvf = this->thermo().p();

    reactionEvaluationScope scope(*this);

    const label nZoneCells = zone_.nCells();
    for(label zci = 0; zci<nZoneCells; zci++)
    {
        const label celli = zone_.celli(zci);

        const scalar rho = rhovf[celli];
        const scalar T = Tvf[celli];
        const scalar p = pvf[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            s_.c[i] = rho*Yvf_[i][celli]/specieWs_[i];
        }

        s_.reset(p, T, celli);

        RR[celli] = reactions_[reactioni].omega(s_);
    }

    return tRR;
}


Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::fvMesh>>
Foam::chemistryModels::standard::specieReactionRR
(
    const label reactioni
) const
{
    PtrList<volInternalScalarField> RR(nSpecie_);
    for (label i=0; i<nSpecie_; i++)
    {
        RR.set
        (
            i,
            volInternalScalarField::New
            (
                "RR:" + reactions_[reactioni].name() + ":" + Yvf_[i].name(),
                this->mesh(),
                dimensionedScalar(dimensions::density/dimensions::time, 0)
            ).ptr()
        );
    }

    if (!this->chemistry_ || reduction_.reactionDisabled(reactioni))
    {
        return RR;
    }

    tmp<volScalarField> trhovf(this->thermo().rho());
    const volScalarField& rhovf = trhovf();

    const volScalarField& Tvf = this->thermo().T();
    const volScalarField& pvf = this->thermo().p();

    scalarField& dNdtByV = YTWork_[0];

    reactionEvaluationScope scope(*this);

    const label nZoneCells = zone_.nCells();
    for(label zci = 0; zci<nZoneCells; zci++)
    {
        const label celli = zone_.celli(zci);

        const scalar rho = rhovf[celli];
        const scalar T = Tvf[celli];
        const scalar p = pvf[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            s_.c[i] = rho*Yvf_[i][celli]/specieWs_[i];
        }

        s_.reset(p, T, celli);

        dNdtByV = Zero;

        reactions_[reactioni].dNdtByV(s_, dNdtByV, reduced_, c2s_);

        for (label i=0; i<nSpecie_; i++)
        {
            RR[i][celli] = dNdtByV[i]*specieWs_[i];
        }
    }

    return RR;
}


void Foam::chemistryModels::standard::calculate()
{
    if (!this->chemistry_)
    {
        return;
    }

    if (!zone_.all())
    {
        forAll(RR_, fieldi)
        {
            RR_[fieldi] = Zero;
        }
    }

    tmp<volScalarField> trhovf(this->thermo().rho());
    const volScalarField& rhovf = trhovf();

    const volScalarField& Tvf = this->thermo().T();
    const volScalarField& pvf = this->thermo().p();

    scalarField& dNdtByV = YTWork_[0];

    reactionEvaluationScope scope(*this);

    const label nZoneCells = zone_.nCells();
    for(label zci = 0; zci<nZoneCells; zci++)
    {
        const label celli = zone_.celli(zci);

        const scalar rho = rhovf[celli];
        const scalar T = Tvf[celli];
        const scalar p = pvf[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            s_.c[i] = rho*Yvf_[i][celli]/specieWs_[i];
        }

        s_.reset(p, T, celli);

        dNdtByV = Zero;

        forAll(reactions_, ri)
        {
            if (!reduction_.reactionDisabled(ri))
            {
                reactions_[ri].dNdtByV(s_, dNdtByV, reduced_, c2s_);
            }
        }

        for (label i=0; i<reduction_.nActiveSpecies(); i++)
        {
            RR_[sToc(i)][celli] = dNdtByV[i]*specieWs_[sToc(i)];
        }
    }
}


Foam::scalar Foam::chemistryModels::standard::solve
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


Foam::scalar Foam::chemistryModels::standard::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


Foam::tmp<Foam::volInternalScalarField>
Foam::chemistryModels::standard::tc() const
{
    tmp<volInternalScalarField> ttc
    (
        volInternalScalarField::New
        (
            "tc",
            this->mesh(),
            dimensionedScalar(dimensions::time, small)
        )
    );
    scalarField& tc = ttc.ref();

    if (!this->chemistry_)
    {
        return ttc;
    }

    tmp<volScalarField> trhovf(this->thermo().rho());
    const volScalarField& rhovf = trhovf();

    const volScalarField& Tvf = this->thermo().T();
    const volScalarField& pvf = this->thermo().p();

    reactionEvaluationScope scope(*this);

    const label nZoneCells = zone_.nCells();
    for(label zci = 0; zci<nZoneCells; zci++)
    {
        const label celli = zone_.celli(zci);

        const scalar rho = rhovf[celli];
        const scalar T = Tvf[celli];
        const scalar p = pvf[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            s_.c[i] = rho*Yvf_[i][celli]/specieWs_[i];
        }

        s_.reset(p, T, celli);

        // A reaction's rate scale is calculated as its molar
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
            scalar omegaf, omegar;
            reactions_[i].omega(s_, omegaf, omegar);

            scalar wf = 0;
            forAll(reactions_[i].rhs(), j)
            {
                wf += reactions_[i].rhs()[j].stoichCoeff*omegaf;
            }
            sumW += wf;
            sumWRateByCTot += sqr(wf);

            scalar wr = 0;
            forAll(reactions_[i].lhs(), j)
            {
                wr += reactions_[i].lhs()[j].stoichCoeff*omegar;
            }
            sumW += wr;
            sumWRateByCTot += sqr(wr);
        }

        tc[celli] =
            sumWRateByCTot == 0 ? vGreat : sumW/sumWRateByCTot*sum(s_.c);
    }

    return ttc;
}


Foam::tmp<Foam::volInternalScalarField>
Foam::chemistryModels::standard::Qdot() const
{
    tmp<volInternalScalarField> tQdot
    (
        volInternalScalarField::New
        (
            "Qdot",
            this->mesh_,
            dimensionedScalar(dimensions::powerDensity, 0)
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
        const label nZoneCells = zone_.nCells();

        for(label zci = 0; zci<nZoneCells; zci++)
        {
            const label celli = zone_.celli(zci);

            Qdot[celli] -= specieHfs_[i]*RR_[i][celli];
        }
    }

    return tQdot;
}


void Foam::chemistryModels::standard::solve
(
    scalar& p,
    scalar& T,
    scalarField& c,
    const label li,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    // Reset the size of the ODE system to the simplified size when mechanism
    // reduction is active
    if (odeSolver_->resize())
    {
        odeSolver_->resizeField(cT_);
    }

    const label nSpecie = this->nSpecie();

    // Copy the concentration, T and P to the total solve-vector
    for (int i=0; i<nSpecie; i++)
    {
        cT_[i] = c[i];
    }
    cT_[nSpecie] = T;

    if (debug)
    {
        scalarField dcT(this->nEqns(), rootSmall);
        dcT[nSpecie] = T*rootSmall;
        this->check(0, cT_, dcT, li);
    }

    odeSolver_->solve(0, deltaT, cT_, li, subDeltaT);

    for (int i=0; i<nSpecie; i++)
    {
        c[i] = max(0.0, cT_[i]);
    }
    T = cT_[nSpecie];
}


// ************************************************************************* //
