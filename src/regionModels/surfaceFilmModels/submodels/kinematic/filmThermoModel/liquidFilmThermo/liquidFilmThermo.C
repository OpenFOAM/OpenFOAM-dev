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

#include "liquidFilmThermo.H"
#include "demandDrivenData.H"
#include "thermoSingleLayer.H"
#include "SLGThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(liquidFilmThermo, 0);

addToRunTimeSelectionTable
(
    filmThermoModel,
    liquidFilmThermo,
    dictionary
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const thermoSingleLayer& liquidFilmThermo::thermoFilm() const
{
    if (!isA<thermoSingleLayer>(owner_))
    {
        FatalErrorIn
        (
            "const thermoSingleLayer& liquidFilmThermo::thermoFilm() const"
        )
            << "Thermo model requires a " << thermoSingleLayer::typeName
            << " film to supply the pressure and temperature, but "
            << owner_.type() << " film model selected.  "
            << "Use the 'useReferenceValues' flag to employ reference "
            << "pressure and temperature" << exit(FatalError);
    }

    return refCast<const thermoSingleLayer>(owner_);
}


void liquidFilmThermo::initLiquid(const dictionary& dict)
{
    if (liquidPtr_ != NULL)
    {
        return;
    }

    dict.lookup("liquid") >> name_;

    if (owner_.primaryMesh().foundObject<SLGThermo>("SLGThermo"))
    {
        // retrieve from film thermo
        ownLiquid_ = false;

        const SLGThermo& thermo =
            owner_.primaryMesh().lookupObject<SLGThermo>("SLGThermo");
        label id = thermo.liquidId(name_);
        liquidPtr_ = &thermo.liquids().properties()[id];
    }
    else
    {
        // new liquid create
        ownLiquid_ = true;

        liquidPtr_ = new liquidProperties(dict.subDict(name_ + "Coeffs"));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

liquidFilmThermo::liquidFilmThermo
(
    surfaceFilmModel& owner,
    const dictionary& dict
)
:
    filmThermoModel(typeName, owner, dict),
    name_("unknown_liquid"),
    liquidPtr_(NULL),
    ownLiquid_(false),
    useReferenceValues_(readBool(coeffDict_.lookup("useReferenceValues"))),
    pRef_(0.0),
    TRef_(0.0)
{
    initLiquid(coeffDict_);

    if (useReferenceValues_)
    {
        coeffDict_.lookup("pRef") >> pRef_;
        coeffDict_.lookup("TRef") >> TRef_;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

liquidFilmThermo::~liquidFilmThermo()
{
    if (ownLiquid_)
    {
        deleteDemandDrivenData(liquidPtr_);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const word& liquidFilmThermo::name() const
{
    return name_;
}


scalar liquidFilmThermo::rho
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->rho(p, T);
}


scalar liquidFilmThermo::mu
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->mu(p, T);
}


scalar liquidFilmThermo::sigma
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->sigma(p, T);
}


scalar liquidFilmThermo::Cp
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->Cp(p, T);
}


scalar liquidFilmThermo::kappa
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->K(p, T);
}


scalar liquidFilmThermo::D
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->D(p, T);
}


scalar liquidFilmThermo::hl
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->hl(p, T);
}


scalar liquidFilmThermo::pv
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->pv(p, T);
}


scalar liquidFilmThermo::W() const
{
    return liquidPtr_->W();
}


scalar liquidFilmThermo::Tb(const scalar p) const
{
    return liquidPtr_->pvInvert(p);
}


tmp<volScalarField> liquidFilmThermo::rho() const
{
    tmp<volScalarField> trho
    (
        new volScalarField
        (
            IOobject
            (
                type() + ":rho",
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

    scalarField& rho = trho().internalField();

    if (useReferenceValues_)
    {
        forAll(rho, cellI)
        {
            rho[cellI] = this->rho(pRef_, TRef_);
        }
    }
    else
    {
        const thermoSingleLayer& film = thermoFilm();

        const volScalarField& T = film.T();
        const volScalarField& p = film.pPrimary();

        forAll(rho, cellI)
        {
            rho[cellI] = this->rho(p[cellI], T[cellI]);
        }
    }

    trho().correctBoundaryConditions();

    return trho;
}


tmp<volScalarField> liquidFilmThermo::mu() const
{
    tmp<volScalarField> tmu
    (
        new volScalarField
        (
            IOobject
            (
                type() + ":mu",
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

    scalarField& mu = tmu().internalField();

    if (useReferenceValues_)
    {
        forAll(mu, cellI)
        {
            mu[cellI] = this->mu(pRef_, TRef_);
        }
    }
    else
    {
        const thermoSingleLayer& film = thermoFilm();

        const volScalarField& T = film.T();
        const volScalarField& p = film.pPrimary();

        forAll(mu, cellI)
        {
            mu[cellI] = this->mu(p[cellI], T[cellI]);
        }
    }

    tmu().correctBoundaryConditions();

    return tmu;
}


tmp<volScalarField> liquidFilmThermo::sigma() const
{
    tmp<volScalarField> tsigma
    (
        new volScalarField
        (
            IOobject
            (
                type() + ":sigma",
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

    scalarField& sigma = tsigma().internalField();

    if (useReferenceValues_)
    {
        forAll(sigma, cellI)
        {
            sigma[cellI] = this->sigma(pRef_, TRef_);
        }
    }
    else
    {
        const thermoSingleLayer& film = thermoFilm();

        const volScalarField& T = film.T();
        const volScalarField& p = film.pPrimary();

        forAll(sigma, cellI)
        {
            sigma[cellI] = this->sigma(p[cellI], T[cellI]);
        }
    }

    tsigma().correctBoundaryConditions();

    return tsigma;
}


tmp<volScalarField> liquidFilmThermo::Cp() const
{
    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                type() + ":Cp",
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

    scalarField& Cp = tCp().internalField();

    if (useReferenceValues_)
    {
        forAll(Cp, cellI)
        {
            Cp[cellI] = this->Cp(pRef_, TRef_);
        }
    }
    else
    {
        const thermoSingleLayer& film = thermoFilm();

        const volScalarField& T = film.T();
        const volScalarField& p = film.pPrimary();

        forAll(Cp, cellI)
        {
            Cp[cellI] = this->Cp(p[cellI], T[cellI]);
        }
    }

    tCp().correctBoundaryConditions();

    return tCp;
}


tmp<volScalarField> liquidFilmThermo::kappa() const
{
    tmp<volScalarField> tkappa
    (
        new volScalarField
        (
            IOobject
            (
                type() + ":kappa",
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

    scalarField& kappa = tkappa().internalField();

    if (useReferenceValues_)
    {
        forAll(kappa, cellI)
        {
            kappa[cellI] = this->kappa(pRef_, TRef_);
        }
    }
    else
    {
        const thermoSingleLayer& film = thermoFilm();

        const volScalarField& T = film.T();
        const volScalarField& p = film.pPrimary();

        forAll(kappa, cellI)
        {
            kappa[cellI] = this->kappa(p[cellI], T[cellI]);
        }
    }

    tkappa().correctBoundaryConditions();

    return tkappa;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
