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

#include "chemistryModel.H"
#include "reactingMixture.H"
#include "UniformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryModel<CompType, ThermoType>::chemistryModel
(
    const fvMesh& mesh
)
:
    CompType(mesh),
    ODESystem(),
    Y_(this->thermo().composition().Y()),
    reactions_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>(this->thermo())
    ),
    specieThermo_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>
            (this->thermo()).speciesData()
    ),

    nSpecie_(Y_.size()),
    nReaction_(reactions_.size()),

    RR_(nSpecie_)
{
    // create the fields for the chemistry sources
    forAll(RR_, fieldI)
    {
        RR_.set
        (
            fieldI,
            new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    "RR." + Y_[fieldI].name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
            )
        );
    }

    Info<< "chemistryModel: Number of species = " << nSpecie_
        << " and reactions = " << nReaction_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryModel<CompType, ThermoType>::~chemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::chemistryModel<CompType, ThermoType>::omega
(
    const scalarField& c,
    const scalar T,
    const scalar p
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    tmp<scalarField> tom(new scalarField(nEqns(), 0.0));
    scalarField& om = tom();

    forAll(reactions_, i)
    {
        const Reaction<ThermoType>& R = reactions_[i];

        scalar omegai = omega
        (
            R, c, T, p, pf, cf, lRef, pr, cr, rRef
        );

        forAll(R.lhs(), s)
        {
            const label si = R.lhs()[s].index;
            const scalar sl = R.lhs()[s].stoichCoeff;
            om[si] -= sl*omegai;
        }

        forAll(R.rhs(), s)
        {
            const label si = R.rhs()[s].index;
            const scalar sr = R.rhs()[s].stoichCoeff;
            om[si] += sr*omegai;
        }
    }

    return tom;
}


template<class CompType, class ThermoType>
Foam::scalar Foam::chemistryModel<CompType, ThermoType>::omegaI
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

    const Reaction<ThermoType>& R = reactions_[index];
    scalar w = omega(R, c, T, p, pf, cf, lRef, pr, cr, rRef);
    return(w);
}


template<class CompType, class ThermoType>
Foam::scalar Foam::chemistryModel<CompType, ThermoType>::omega
(
    const Reaction<ThermoType>& R,
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
    scalarField c2(nSpecie_, 0.0);
    for (label i = 0; i < nSpecie_; i++)
    {
        c2[i] = max(0.0, c[i]);
    }

    const scalar kf = R.kf(p, T, c2);
    const scalar kr = R.kr(kf, p, T, c2);

    pf = 1.0;
    pr = 1.0;

    const label Nl = R.lhs().size();
    const label Nr = R.rhs().size();

    label slRef = 0;
    lRef = R.lhs()[slRef].index;

    pf = kf;
    for (label s = 1; s < Nl; s++)
    {
        const label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            const scalar exp = R.lhs()[slRef].exponent;
            pf *= pow(max(0.0, c[lRef]), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            const scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(0.0, c[si]), exp);
        }
    }
    cf = max(0.0, c[lRef]);

    {
        const scalar exp = R.lhs()[slRef].exponent;
        if (exp < 1.0)
        {
            if (cf > SMALL)
            {
                pf *= pow(cf, exp - 1.0);
            }
            else
            {
                pf = 0.0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1.0);
        }
    }

    label srRef = 0;
    rRef = R.rhs()[srRef].index;

    // find the matrix element and element position for the rhs
    pr = kr;
    for (label s = 1; s < Nr; s++)
    {
        const label si = R.rhs()[s].index;
        if (c[si] < c[rRef])
        {
            const scalar exp = R.rhs()[srRef].exponent;
            pr *= pow(max(0.0, c[rRef]), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            const scalar exp = R.rhs()[s].exponent;
            pr *= pow(max(0.0, c[si]), exp);
        }
    }
    cr = max(0.0, c[rRef]);

    {
        const scalar exp = R.rhs()[srRef].exponent;
        if (exp < 1.0)
        {
            if (cr>SMALL)
            {
                pr *= pow(cr, exp - 1.0);
            }
            else
            {
                pr = 0.0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1.0);
        }
    }

    return pf*cf - pr*cr;
}


template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::derivatives
(
    const scalar time,
    const scalarField &c,
    scalarField& dcdt
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    dcdt = omega(c, T, p);

    // constant pressure
    // dT/dt = ...
    scalar rho = 0.0;
    scalar cSum = 0.0;
    for (label i = 0; i < nSpecie_; i++)
    {
        const scalar W = specieThermo_[i].W();
        cSum += c[i];
        rho += W*c[i];
    }
    scalar cp = 0.0;
    for (label i=0; i<nSpecie_; i++)
    {
        cp += c[i]*specieThermo_[i].cp(p, T);
    }
    cp /= rho;

    scalar dT = 0.0;
    for (label i = 0; i < nSpecie_; i++)
    {
        const scalar hi = specieThermo_[i].ha(p, T);
        dT += hi*dcdt[i];
    }
    dT /= rho*cp;

    dcdt[nSpecie_] = -dT;

    // dp/dt = ...
    dcdt[nSpecie_ + 1] = 0.0;
}


template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::jacobian
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
    forAll(c2, i)
    {
        c2[i] = max(c[i], 0.0);
    }

    for (label i=0; i<nEqns(); i++)
    {
        for (label j=0; j<nEqns(); j++)
        {
            dfdc[i][j] = 0.0;
        }
    }

    // length of the first argument must be nSpecie()
    dcdt = omega(c2, T, p);

    forAll(reactions_, ri)
    {
        const Reaction<ThermoType>& R = reactions_[ri];

        const scalar kf0 = R.kf(p, T, c2);
        const scalar kr0 = R.kr(p, T, c2);

        forAll(R.lhs(), j)
        {
            const label sj = R.lhs()[j].index;
            scalar kf = kf0;
            forAll(R.lhs(), i)
            {
                const label si = R.lhs()[i].index;
                const scalar el = R.lhs()[i].exponent;
                if (i == j)
                {
                    if (el < 1.0)
                    {
                        if (c2[si] > SMALL)
                        {
                            kf *= el*pow(c2[si] + VSMALL, el - 1.0);
                        }
                        else
                        {
                            kf = 0.0;
                        }
                    }
                    else
                    {
                        kf *= el*pow(c2[si], el - 1.0);
                    }
                }
                else
                {
                    kf *= pow(c2[si], el);
                }
            }

            forAll(R.lhs(), i)
            {
                const label si = R.lhs()[i].index;
                const scalar sl = R.lhs()[i].stoichCoeff;
                dfdc[si][sj] -= sl*kf;
            }
            forAll(R.rhs(), i)
            {
                const label si = R.rhs()[i].index;
                const scalar sr = R.rhs()[i].stoichCoeff;
                dfdc[si][sj] += sr*kf;
            }
        }

        forAll(R.rhs(), j)
        {
            const label sj = R.rhs()[j].index;
            scalar kr = kr0;
            forAll(R.rhs(), i)
            {
                const label si = R.rhs()[i].index;
                const scalar er = R.rhs()[i].exponent;
                if (i == j)
                {
                    if (er < 1.0)
                    {
                        if (c2[si] > SMALL)
                        {
                            kr *= er*pow(c2[si] + VSMALL, er - 1.0);
                        }
                        else
                        {
                            kr = 0.0;
                        }
                    }
                    else
                    {
                        kr *= er*pow(c2[si], er - 1.0);
                    }
                }
                else
                {
                    kr *= pow(c2[si], er);
                }
            }

            forAll(R.lhs(), i)
            {
                const label si = R.lhs()[i].index;
                const scalar sl = R.lhs()[i].stoichCoeff;
                dfdc[si][sj] += sl*kr;
            }
            forAll(R.rhs(), i)
            {
                const label si = R.rhs()[i].index;
                const scalar sr = R.rhs()[i].stoichCoeff;
                dfdc[si][sj] -= sr*kr;
            }
        }
    }

    // Calculate the dcdT elements numerically
    const scalar delta = 1.0e-3;
    const scalarField dcdT0(omega(c2, T - delta, p));
    const scalarField dcdT1(omega(c2, T + delta, p));

    for (label i = 0; i < nEqns(); i++)
    {
        dfdc[i][nSpecie()] = 0.5*(dcdT1[i] - dcdT0[i])/delta;
    }
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::chemistryModel<CompType, ThermoType>::tc() const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

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
        this->thermo().rho()
    );

    tmp<volScalarField> ttc
    (
        new volScalarField
        (
            IOobject
            (
                "tc",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimTime, SMALL),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    scalarField& tc = ttc();
    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    const label nReaction = reactions_.size();

    if (this->chemistry_)
    {
        forAll(rho, celli)
        {
            scalar rhoi = rho[celli];
            scalar Ti = T[celli];
            scalar pi = p[celli];
            scalarField c(nSpecie_);
            scalar cSum = 0.0;

            for (label i=0; i<nSpecie_; i++)
            {
                scalar Yi = Y_[i][celli];
                c[i] = rhoi*Yi/specieThermo_[i].W();
                cSum += c[i];
            }

            forAll(reactions_, i)
            {
                const Reaction<ThermoType>& R = reactions_[i];

                omega(R, c, Ti, pi, pf, cf, lRef, pr, cr, rRef);

                forAll(R.rhs(), s)
                {
                    scalar sr = R.rhs()[s].stoichCoeff;
                    tc[celli] += sr*pf*cf;
                }
            }
            tc[celli] = nReaction*cSum/tc[celli];
        }
    }


    ttc().correctBoundaryConditions();

    return ttc;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::chemistryModel<CompType, ThermoType>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );


    if (this->chemistry_)
    {
        scalarField& Sh = tSh();

        forAll(Y_, i)
        {
            forAll(Sh, cellI)
            {
                const scalar hi = specieThermo_[i].Hc();
                Sh[cellI] -= hi*RR_[i][cellI];
            }
        }
    }

    return tSh;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::chemistryModel<CompType, ThermoType>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("dQ", dimEnergy/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        volScalarField& dQ = tdQ();
        dQ.dimensionedInternalField() = this->mesh_.V()*Sh()();
    }

    return tdQ;
}


template<class CompType, class ThermoType>
Foam::label Foam::chemistryModel<CompType, ThermoType>::nEqns() const
{
    // nEqns = number of species + temperature + pressure
    return nSpecie_ + 2;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::chemistryModel<CompType, ThermoType>::calculateRR
(
    const label reactionI,
    const label specieI
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

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
        this->thermo().rho()
    );

    tmp<DimensionedField<scalar, volMesh> > tRR
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "RR",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    DimensionedField<scalar, volMesh>& RR = tRR();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T[celli];
        const scalar pi = p[celli];

        scalarField c(nSpecie_, 0.0);
        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Y_[i][celli];
            c[i] = rhoi*Yi/specieThermo_[i].W();
        }

        const scalar w = omegaI
        (
            reactionI,
            c,
            Ti,
            pi,
            pf,
            cf,
            lRef,
            pr,
            cr,
            rRef
        );

        RR[celli] = w*specieThermo_[specieI].W();

    }

    return tRR;
}


template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::calculate()
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
        this->thermo().rho()
    );

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T[celli];
        const scalar pi = p[celli];

        scalarField c(nSpecie_, 0.0);
        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Y_[i][celli];
            c[i] = rhoi*Yi/specieThermo_[i].W();
        }

        const scalarField dcdt(omega(c, Ti, pi));

        for (label i=0; i<nSpecie_; i++)
        {
            RR_[i][celli] = dcdt[i]*specieThermo_[i].W();
        }
    }
}


template<class CompType, class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::chemistryModel<CompType, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    CompType::correct();

    scalar deltaTMin = GREAT;

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
        this->thermo().rho()
    );

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    scalarField c(nSpecie_);
    scalarField c0(nSpecie_);

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        scalar pi = p[celli];
        scalar Ti = T[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            c[i] = rhoi*Y_[i][celli]/specieThermo_[i].W();
            c0[i] = c[i];
        }

        // Initialise time progress
        scalar timeLeft = deltaT[celli];

        // Calculate the chemical source terms
        while (timeLeft > SMALL)
        {
            scalar dt = timeLeft;
            this->solve(c, Ti, pi, dt, this->deltaTChem_[celli]);
            timeLeft -= dt;
        }

        deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

        for (label i=0; i<nSpecie_; i++)
        {
            RR_[i][celli] = (c[i] - c0[i])*specieThermo_[i].W()/deltaT[celli];
        }
    }

    return deltaTMin;
}


template<class CompType, class ThermoType>
Foam::scalar Foam::chemistryModel<CompType, ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar> >(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}


template<class CompType, class ThermoType>
Foam::scalar Foam::chemistryModel<CompType, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::solve
(
    scalarField &c,
    scalar& T,
    scalar& p,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    notImplemented
    (
        "chemistryModel::solve"
        "("
            "scalarField&, "
            "scalar&, "
            "scalar&, "
            "scalar&, "
            "scalar&"
        ") const"
    );
}


// ************************************************************************* //
