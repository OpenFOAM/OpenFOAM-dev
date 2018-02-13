/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "pyrolysisChemistryModel.H"
#include "solidReaction.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class SolidThermo, class GasThermo>
Foam::pyrolysisChemistryModel<CompType, SolidThermo, GasThermo>::
pyrolysisChemistryModel
(
    typename CompType::reactionThermo& thermo
)
:
    solidChemistryModel<CompType, SolidThermo>(thermo),
    pyrolisisGases_(this->reactions_[0].gasSpecies()),
    gasThermo_(pyrolisisGases_.size()),
    nGases_(pyrolisisGases_.size()),
    nSpecie_(this->Ys_.size() + nGases_),
    RRg_(nGases_),
    Ys0_(this->nSolids_),
    cellCounter_(0)
{
    // create the fields for the chemistry sources
    forAll(this->RRs_, fieldi)
    {
        IOobject header
        (
            this->Ys_[fieldi].name() + "0",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ
        );

        // check if field exists and can be read
        if (header.typeHeaderOk<volScalarField>(true))
        {
            Ys0_.set
            (
                fieldi,
                new volScalarField
                (
                    IOobject
                    (
                        this->Ys_[fieldi].name() + "0",
                        this->mesh().time().timeName(),
                        this->mesh(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    this->mesh()
                )
            );
        }
        else
        {
            volScalarField Y0Default
            (
                IOobject
                (
                    "Y0Default",
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh()
            );

            Ys0_.set
            (
                fieldi,
                new volScalarField
                (
                    IOobject
                    (
                        this->Ys_[fieldi].name() + "0",
                        this->mesh().time().timeName(),
                        this->mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    Y0Default
                )
            );

            // Calculate initial values of Ysi0 = rho*delta*Yi
            Ys0_[fieldi].primitiveFieldRef() =
                this->solidThermo().rho()
               *max(this->Ys_[fieldi], scalar(0.001))*this->mesh().V();
        }
    }

    forAll(RRg_, fieldi)
    {
        RRg_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RRg." + pyrolisisGases_[fieldi],
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
            )
        );
    }

    forAll(gasThermo_, gasI)
    {
        dictionary thermoDict =
            this->mesh().template lookupObject<dictionary>
            (
                basicThermo::dictName
            ).subDict(pyrolisisGases_[gasI]);

        gasThermo_.set
        (
            gasI,
            new GasThermo(thermoDict)
        );
    }

    Info<< "pyrolysisChemistryModel: " << nl;
    Info<< indent << "Number of solids = " << this->nSolids_ << nl;
    Info<< indent << "Number of gases = " << nGases_ << nl;
    forAll(this->reactions_, i)
    {
        Info<< dynamic_cast<const solidReaction<SolidThermo>& >
        (
            this->reactions_[i]
        ) << nl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class SolidThermo, class GasThermo>
Foam::pyrolysisChemistryModel<CompType, SolidThermo, GasThermo>::
~pyrolysisChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class SolidThermo, class GasThermo>
Foam::scalarField Foam::
pyrolysisChemistryModel<CompType, SolidThermo, GasThermo>::omega
(
    const scalarField& c,
    const scalar T,
    const scalar p,
    const bool updateC0
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    const label celli = cellCounter_;

    scalarField om(nEqns(), 0.0);

    forAll(this->reactions_, i)
    {
        const Reaction<SolidThermo>& R = this->reactions_[i];

        scalar omegai = omega
        (
            R, c, T, p, pf, cf, lRef, pr, cr, rRef
        );
        scalar rhoL = 0.0;
        forAll(R.lhs(), s)
        {
            label si = R.lhs()[s].index;
            om[si] -= omegai;
            rhoL = this->solidThermo_[si].rho(p, T);
        }
        scalar sr = 0.0;
        forAll(R.rhs(), s)
        {
            label si = R.rhs()[s].index;
            scalar rhoR = this->solidThermo_[si].rho(p, T);
            sr = rhoR/rhoL;
            om[si] += sr*omegai;

            if (updateC0)
            {
                Ys0_[si][celli] += sr*omegai;
            }
        }
        forAll(R.grhs(), g)
        {
            label gi = R.grhs()[g].index;
            om[gi + this->nSolids_] += (1.0 - sr)*omegai;
        }
    }

    return om;
}


template<class CompType, class SolidThermo, class GasThermo>
Foam::scalar
Foam::pyrolysisChemistryModel<CompType, SolidThermo, GasThermo>::omega
(
    const Reaction<SolidThermo>& R,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    scalarField c1(nSpecie_, 0.0);

    label celli = cellCounter_;

    for (label i=0; i<nSpecie_; i++)
    {
        c1[i] = max(0.0, c[i]);
    }

    scalar kf = R.kf(p, T, c1);

    const label Nl = R.lhs().size();

    for (label s=0; s<Nl; s++)
    {
        label si = R.lhs()[s].index;
        const scalar exp = R.lhs()[si].exponent;

        kf *=
            pow(c1[si]/Ys0_[si][celli], exp)
           *(Ys0_[si][celli]);
    }

    return kf;
}


template<class CompType, class SolidThermo, class GasThermo>
Foam::scalar Foam::pyrolysisChemistryModel<CompType, SolidThermo, GasThermo>::
omegaI
(
    const label index,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{

    const Reaction<SolidThermo>& R = this->reactions_[index];
    scalar w = omega(R, c, T, p, pf, cf, lRef, pr, cr, rRef);
    return(w);
}


template<class CompType, class SolidThermo, class GasThermo>
void Foam::pyrolysisChemistryModel<CompType, SolidThermo, GasThermo>::
derivatives
(
    const scalar time,
    const scalarField &c,
    scalarField& dcdt
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    dcdt = 0.0;

    dcdt = omega(c, T, p);

    //Total mass concentration
    scalar cTot = 0.0;
    for (label i=0; i<this->nSolids_; i++)
    {
        cTot += c[i];
    }

    scalar newCp = 0.0;
    scalar newhi = 0.0;
    for (label i=0; i<this->nSolids_; i++)
    {
        scalar dYidt = dcdt[i]/cTot;
        scalar Yi = c[i]/cTot;
        newCp += Yi*this->solidThermo_[i].Cp(p, T);
        newhi -= dYidt*this->solidThermo_[i].Hc();
    }

    scalar dTdt = newhi/newCp;
    scalar dtMag = min(500.0, mag(dTdt));
    dcdt[nSpecie_] = dTdt*dtMag/(mag(dTdt) + 1.0e-10);

    // dp/dt = ...
    dcdt[nSpecie_ + 1] = 0.0;
}


template<class CompType, class SolidThermo, class GasThermo>
void Foam::pyrolysisChemistryModel<CompType, SolidThermo, GasThermo>::
jacobian
(
    const scalar t,
    const scalarField& c,
    scalarField& dcdt,
    scalarSquareMatrix& dfdc
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    scalarField c2(nSpecie_, 0.0);

    for (label i=0; i<this->nSolids_; i++)
    {
        c2[i] = max(c[i], 0.0);
    }

    for (label i=0; i<nEqns(); i++)
    {
        for (label j=0; j<nEqns(); j++)
        {
            dfdc(i, j) = 0.0;
        }
    }

    // length of the first argument must be nSolids
    dcdt = omega(c2, T, p);

    for (label ri=0; ri<this->reactions_.size(); ri++)
    {
        const Reaction<SolidThermo>& R = this->reactions_[ri];

        scalar kf0 = R.kf(p, T, c2);

        forAll(R.lhs(), j)
        {
            label sj = R.lhs()[j].index;
            scalar kf = kf0;
            forAll(R.lhs(), i)
            {
                label si = R.lhs()[i].index;
                scalar exp = R.lhs()[i].exponent;
                if (i == j)
                {
                    if (exp < 1.0)
                    {
                        if (c2[si] > small)
                        {
                            kf *= exp*pow(c2[si], exp - 1.0);
                        }
                        else
                        {
                            kf = 0.0;
                        }
                    }
                    else
                    {
                        kf *= exp*pow(c2[si], exp - 1.0);
                    }
                }
                else
                {
                    Info<< "Solid reactions have only elements on slhs"
                        << endl;
                    kf = 0.0;
                }
            }

            forAll(R.lhs(), i)
            {
                label si = R.lhs()[i].index;
                dfdc[si][sj] -= kf;
            }
            forAll(R.rhs(), i)
            {
                label si = R.rhs()[i].index;
                dfdc[si][sj] += kf;
            }
        }
    }

    // calculate the dcdT elements numerically
    scalar delta = 1.0e-8;
    scalarField dcdT0 = omega(c2, T - delta, p);
    scalarField dcdT1 = omega(c2, T + delta, p);

    for (label i=0; i<nEqns(); i++)
    {
        dfdc[i][nSpecie_] = 0.5*(dcdT1[i] - dcdT0[i])/delta;
    }

}


template<class CompType, class SolidThermo, class GasThermo>
Foam::label Foam::
pyrolysisChemistryModel<CompType, SolidThermo, GasThermo>::nEqns() const
{
    // nEqns = number of solids + gases + temperature + pressure
    return (nSpecie_ + 2);
}


template<class CompType, class SolidThermo, class GasThermo>
void Foam::pyrolysisChemistryModel<CompType, SolidThermo, GasThermo>::
calculate()
{
    if (!this->chemistry_)
    {
        return;
    }

    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->solidThermo().rho()
    );

    forAll(this->RRs_, i)
    {
        this->RRs_[i].field() = 0.0;
    }

    forAll(RRg_, i)
    {
        RRg_[i].field() = 0.0;
    }

    forAll(rho, celli)
    {
        cellCounter_ = celli;

        const scalar delta = this->mesh().V()[celli];

        if (this->reactingCells_[celli])
        {
            scalar rhoi = rho[celli];
            scalar Ti = this->solidThermo().T()[celli];
            scalar pi = this->solidThermo().p()[celli];

            scalarField c(nSpecie_, 0.0);
            for (label i=0; i<this->nSolids_; i++)
            {
                c[i] = rhoi*this->Ys_[i][celli]*delta;
            }

            const scalarField dcdt = omega(c, Ti, pi, true);

            forAll(this->RRs_, i)
            {
                this->RRs_[i][celli] = dcdt[i]/delta;
            }

            forAll(RRg_, i)
            {
                RRg_[i][celli] = dcdt[this->nSolids_ + i]/delta;
            }
        }
    }
}


template<class CompType, class SolidThermo, class GasThermo>
Foam::scalar
Foam::pyrolysisChemistryModel<CompType, SolidThermo, GasThermo>::solve
(
    const scalar deltaT
)
{
    scalar deltaTMin = great;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->solidThermo().rho()
    );

    forAll(this->RRs_, i)
    {
        this->RRs_[i].field() = 0.0;
    }
    forAll(RRg_, i)
    {
        RRg_[i].field() = 0.0;
    }

    const scalarField& T = this->solidThermo().T();
    const scalarField& p = this->solidThermo().p();

    scalarField c(nSpecie_, 0.0);
    scalarField c0(nSpecie_, 0.0);
    scalarField dc(nSpecie_, 0.0);
    scalarField delta(this->mesh().V());

    forAll(rho, celli)
    {
        if (this->reactingCells_[celli])
        {
            cellCounter_ = celli;

            scalar rhoi = rho[celli];
            scalar pi = p[celli];
            scalar Ti = T[celli];

            for (label i=0; i<this->nSolids_; i++)
            {
                c[i] = rhoi*this->Ys_[i][celli]*delta[celli];
            }

            c0 = c;

            // Initialise time progress
            scalar timeLeft = deltaT;

            // calculate the chemical source terms
            while (timeLeft > small)
            {
                scalar dt = timeLeft;
                this->solve(c, Ti, pi, dt, this->deltaTChem_[celli]);
                timeLeft -= dt;
            }

            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);
            dc = c - c0;

            forAll(this->RRs_, i)
            {
                this->RRs_[i][celli] = dc[i]/(deltaT*delta[celli]);
            }

            forAll(RRg_, i)
            {
                RRg_[i][celli] = dc[this->nSolids_ + i]/(deltaT*delta[celli]);
            }

            // Update Ys0_
            dc = omega(c0, Ti, pi, true);
        }
    }

    // Don't allow the time-step to change more than a factor of 2
    deltaTMin = min(deltaTMin, 2*deltaT);

    return deltaTMin;
}


template<class CompType, class SolidThermo,class GasThermo>
Foam::tmp<Foam::volScalarField>
Foam::pyrolysisChemistryModel<CompType, SolidThermo, GasThermo>::gasHs
(
    const volScalarField& p,
    const volScalarField& T,
    const label index
) const
{
    tmp<volScalarField> tHs
    (
        new volScalarField
        (
            IOobject
            (
                "Hs_" + pyrolisisGases_[index],
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimMass, 0.0)
        )
    );

    volScalarField::Internal& gasHs = tHs.ref();

    const GasThermo& mixture = gasThermo_[index];

    forAll(gasHs, celli)
    {
        gasHs[celli] = mixture.Hs(p[celli], T[celli]);
    }

    return tHs;
}


template<class CompType, class SolidThermo, class GasThermo>
void Foam::pyrolysisChemistryModel<CompType, SolidThermo, GasThermo>::solve
(
    scalarField &c,
    scalar& T,
    scalar& p,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    NotImplemented;
}


// ************************************************************************* //
