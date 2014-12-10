/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "constantFilmThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(constantFilmThermo, 0);

addToRunTimeSelectionTable
(
    filmThermoModel,
    constantFilmThermo,
    dictionary
);


void constantFilmThermo::init(thermoData& td)
{
    if (coeffDict_.readIfPresent(td.name_, td.value_))
    {
        td.set_ = true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constantFilmThermo::constantFilmThermo
(
    surfaceFilmModel& owner,
    const dictionary& dict
)
:
    filmThermoModel(typeName, owner, dict),
    name_(coeffDict_.lookup("specieName")),
    rho0_("rho0"),
    mu0_("mu0"),
    sigma0_("sigma0"),
    Cp0_("Cp0"),
    kappa0_("kappa0"),
    hl0_("hl0"),
    pv0_("pv0"),
    W0_("W0"),
    Tb0_("Tb0")
{
    init(rho0_);
    init(mu0_);
    init(sigma0_);
    init(Cp0_);
    init(kappa0_);
    init(hl0_);
    init(pv0_);
    init(W0_);
    init(Tb0_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constantFilmThermo::~constantFilmThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const word& constantFilmThermo::name() const
{
    return name_;
}


scalar constantFilmThermo::rho
(
    const scalar p,
    const scalar T
) const
{
    if (!rho0_.set_)
    {
        coeffDict_.lookup(rho0_.name_) >> rho0_.value_;
        rho0_.set_ = true;
    }

    return rho0_.value_;
}


scalar constantFilmThermo::mu
(
    const scalar p,
    const scalar T
) const
{
    if (!mu0_.set_)
    {
        coeffDict_.lookup(mu0_.name_) >> mu0_.value_;
        mu0_.set_ = true;
    }

    return mu0_.value_;
}


scalar constantFilmThermo::sigma
(
    const scalar p,
    const scalar T
) const
{
    if (!sigma0_.set_)
    {
        coeffDict_.lookup(sigma0_.name_) >> sigma0_.value_;
        sigma0_.set_ = true;
    }

    return sigma0_.value_;
}


scalar constantFilmThermo::Cp
(
    const scalar p,
    const scalar T
) const
{
    if (!Cp0_.set_)
    {
        coeffDict_.lookup(Cp0_.name_) >> Cp0_.value_;
        Cp0_.set_ = true;
    }

    return Cp0_.value_;
}


scalar constantFilmThermo::kappa
(
    const scalar p,
    const scalar T
) const
{
    if (!kappa0_.set_)
    {
        coeffDict_.lookup(kappa0_.name_) >> kappa0_.value_;
        kappa0_.set_ = true;
    }

    return kappa0_.value_;
}


scalar constantFilmThermo::D
(
    const scalar p,
    const scalar T
) const
{
    if (!D0_.set_)
    {
        coeffDict_.lookup(D0_.name_) >> D0_.value_;
        D0_.set_ = true;
    }

    return D0_.value_;
}


scalar constantFilmThermo::hl
(
    const scalar p,
    const scalar T
) const
{
    if (!hl0_.set_)
    {
        coeffDict_.lookup(hl0_.name_) >> hl0_.value_;
        hl0_.set_ = true;
    }

    return hl0_.value_;
}


scalar constantFilmThermo::pv
(
    const scalar p,
    const scalar T
) const
{
    if (!pv0_.set_)
    {
        coeffDict_.lookup(pv0_.name_) >> pv0_.value_;
        pv0_.set_ = true;
    }

    return pv0_.value_;
}


scalar constantFilmThermo::W() const
{
    if (!W0_.set_)
    {
        coeffDict_.lookup(W0_.name_) >> W0_.value_;
        W0_.set_ = true;
    }

    return W0_.value_;
}


scalar constantFilmThermo::Tb(const scalar p) const
{
    if (!Tb0_.set_)
    {
        coeffDict_.lookup(Tb0_.name_) >> Tb0_.value_;
        Tb0_.set_ = true;
    }

    return Tb0_.value_;
}


tmp<volScalarField> constantFilmThermo::rho() const
{
    tmp<volScalarField> trho
    (
        new volScalarField
        (
            IOobject
            (
                type() + ':' + rho0_.name_,
                owner().time().timeName(),
                owner().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            owner().regionMesh(),
            dimensionedScalar("0", dimDensity, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    trho().internalField() = this->rho(0, 0);
    trho().correctBoundaryConditions();

    return trho;
}


tmp<volScalarField> constantFilmThermo::mu() const
{
    tmp<volScalarField> tmu
    (
        new volScalarField
        (
            IOobject
            (
                type() + ':' + mu0_.name_,
                owner().time().timeName(),
                owner().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            owner().regionMesh(),
            dimensionedScalar("0", dimPressure*dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tmu().internalField() = this->mu(0, 0);
    tmu().correctBoundaryConditions();

    return tmu;
}


tmp<volScalarField> constantFilmThermo::sigma() const
{
    tmp<volScalarField> tsigma
    (
        new volScalarField
        (
            IOobject
            (
                type() + ':' + sigma0_.name_,
                owner().time().timeName(),
                owner().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            owner().regionMesh(),
            dimensionedScalar("0", dimMass/sqr(dimTime), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tsigma().internalField() = this->sigma(0, 0);
    tsigma().correctBoundaryConditions();

    return tsigma;
}


tmp<volScalarField> constantFilmThermo::Cp() const
{
    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                type() + ':' + Cp0_.name_,
                owner().time().timeName(),
                owner().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            owner().regionMesh(),
            dimensionedScalar("0", dimEnergy/dimMass/dimTemperature, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tCp().internalField() = this->Cp(0, 0);
    tCp().correctBoundaryConditions();

    return tCp;
}


tmp<volScalarField> constantFilmThermo::kappa() const
{
    tmp<volScalarField> tkappa
    (
        new volScalarField
        (
            IOobject
            (
                type() + ':' + kappa0_.name_,
                owner().time().timeName(),
                owner().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            owner().regionMesh(),
            dimensionedScalar("0", dimPower/dimLength/dimTemperature, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tkappa().internalField() = this->kappa(0, 0);
    tkappa().correctBoundaryConditions();

    return tkappa;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
