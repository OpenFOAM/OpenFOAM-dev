/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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

#include "ubRhoThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ubRhoThermo, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ubRhoThermo::ubRhoThermo(const fvMesh& mesh)
:
    unburntPhaseName_("u"),
    burntPhaseName_("b"),
    b_
    (
        IOobject
        (
            "b",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    c_("c", scalar(1) - b_),
    uThermo_(uRhoMulticomponentThermo::New(mesh, unburntPhaseName_)),
    bThermo_(bRhoMulticomponentThermo::New(mesh, burntPhaseName_)),
    rho_("rho", 1.0/(b_/uThermo_->rho() + c_/bThermo_->rho())),
    psi_("psi", 1.0/(b_/uThermo_->psi() + c_/bThermo_->psi())),
    mu_("mu", b_*uThermo_->mu() + c_*bThermo_->mu()),
    kappa_("kappa", b_*uThermo_->kappa() + c_*bThermo_->kappa()),
    alphau_
    (
        phasePropertyName("alpha", unburntPhaseName_), rho_*b_/uThermo_->rho()
    ),
    alphab_
    (
        phasePropertyName("alpha", burntPhaseName_), rho_*c_/bThermo_->rho()
    )
{
    uThermo_->validate
    (
        IOobject::groupName(type(), unburntPhaseName_),
        "h"
    );

    bThermo_->validate
    (
        IOobject::groupName(type(), burntPhaseName_),
        "h"
    );

    b_.oldTime();
    c_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ubRhoThermo::~ubRhoThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::IOdictionary& Foam::ubRhoThermo::properties() const
{
    NotImplemented;
    return uThermo_->properties();
}


Foam::IOdictionary& Foam::ubRhoThermo::properties()
{
    NotImplemented;
    return uThermo_->properties();
}


const Foam::fvMesh& Foam::ubRhoThermo::mesh() const
{
    return uThermo_->mesh();
}


const Foam::word& Foam::ubRhoThermo::phaseName() const
{
    return word::null;
}


Foam::word Foam::ubRhoThermo::mixtureName() const
{
    NotImplemented;
    return word::null;
}


Foam::word Foam::ubRhoThermo::thermoName() const
{
    NotImplemented;
    return word::null;
}


void Foam::ubRhoThermo::correct()
{
    uThermo_->correct();
    bThermo_->correct();

    rho_ = 1.0/(b_/uThermo_->rho() + c_/bThermo_->rho());
    psi_ = 1.0/(b_/uThermo_->psi() + c_/bThermo_->psi());
    mu_ = b_*uThermo_->mu() + c_*bThermo_->mu();
    kappa_ = b_*uThermo_->kappa() + c_*bThermo_->kappa();

    alphau_ = rho_*b_/uThermo_->rho();
    alphab_ = rho_*c_/bThermo_->rho();
}


void Foam::ubRhoThermo::reset()
{
    if (uThermo_->containsSpecie("egr"))
    {
        uThermo_->reset(b_, c_, bThermo_->Y(), bThermo_->he());
    }
    else
    {
        FatalErrorInFunction
            << "EGR not supported by " << uThermo_->type()
            << exit(FatalError);
    }
}


Foam::tmp<Foam::volScalarField> Foam::ubRhoThermo::W() const
{
    NotImplemented;
    return uThermo_->W();
}


Foam::tmp<Foam::scalarField> Foam::ubRhoThermo::W(const label patchi) const
{
    NotImplemented;
    return uThermo_->W(patchi);
}


const Foam::volScalarField& Foam::ubRhoThermo::p() const
{
    return uThermo_->p();
}


Foam::volScalarField& Foam::ubRhoThermo::p()
{
    return uThermo_->p();
}


const Foam::volScalarField& Foam::ubRhoThermo::psi() const
{
    return psi_;
}


const Foam::volScalarField& Foam::ubRhoThermo::T() const
{
    NotImplemented;
    return uThermo_->T();
}


Foam::volScalarField& Foam::ubRhoThermo::T()
{
    NotImplemented;
    return uThermo_->T();
}

const Foam::volScalarField& Foam::ubRhoThermo::he() const
{
    NotImplemented;
    return uThermo_->he();
}

Foam::volScalarField& Foam::ubRhoThermo::he()
{
    NotImplemented;
    return uThermo_->he();
}

const Foam::volScalarField& Foam::ubRhoThermo::Cp() const
{
    NotImplemented;
    return uThermo_->Cp();
}

const Foam::volScalarField& Foam::ubRhoThermo::Cv() const
{
    NotImplemented;
    return uThermo_->Cv();
}

const Foam::volScalarField& Foam::ubRhoThermo::Cpv() const
{
    NotImplemented;
    return uThermo_->Cpv();
}


Foam::tmp<Foam::volScalarField> Foam::ubRhoThermo::rho() const
{
    return rho_;
}


Foam::tmp<Foam::scalarField> Foam::ubRhoThermo::rho(const label patchi) const
{
    return rho_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::ubRhoThermo::rho()
{
    return rho_;
}


Foam::tmp<Foam::volScalarField> Foam::ubRhoThermo::he
(
    const Foam::volScalarField& p,
    const Foam::volScalarField& T
) const
{
    NotImplemented;
    return uThermo_->he(p, T);
}


Foam::tmp<Foam::volScalarField::Internal> Foam::ubRhoThermo::he
(
    const Foam::volScalarField::Internal& p,
    const Foam::volScalarField::Internal& T
) const
{
    NotImplemented;
    return uThermo_->he(p, T);
}


Foam::tmp<Foam::scalarField> Foam::ubRhoThermo::he
(
    const Foam::scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return uThermo_->he(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::ubRhoThermo::he
(
    const Foam::scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return uThermo_->he(T, patchi);
}


Foam::tmp<Foam::volScalarField::Internal> Foam::ubRhoThermo::he
(
    const Foam::volScalarField::Internal& T,
    const fvSource& model,
    const Foam::volScalarField::Internal& source
) const
{
    NotImplemented;
    return uThermo_->he(T, model, source);
}


Foam::tmp<Foam::scalarField> Foam::ubRhoThermo::he
(
    const Foam::scalarField& T,
    const fvSource& model,
    const Foam::scalarField& source,
    const labelUList& cells
) const
{
    NotImplemented;
    return uThermo_->he(T, model, source, cells);
}


Foam::tmp<Foam::volScalarField> Foam::ubRhoThermo::hs() const
{
    NotImplemented;
    return uThermo_->hs();
}

Foam::tmp<Foam::volScalarField> Foam::ubRhoThermo::hs
(
    const Foam::volScalarField& p,
    const Foam::volScalarField& T
) const
{
    NotImplemented;
    return uThermo_->hs(p, T);
}


Foam::tmp<Foam::volScalarField::Internal> Foam::ubRhoThermo::hs
(
    const Foam::volScalarField::Internal& p,
    const Foam::volScalarField::Internal& T
) const
{
    NotImplemented;
    return uThermo_->hs(p, T);
}


Foam::tmp<Foam::scalarField> Foam::ubRhoThermo::hs
(
    const Foam::scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return uThermo_->hs(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::ubRhoThermo::hs
(
    const Foam::scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return uThermo_->hs(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::ubRhoThermo::ha() const
{
    NotImplemented;
    return uThermo_->ha();
}


Foam::tmp<Foam::volScalarField> Foam::ubRhoThermo::ha
(
    const Foam::volScalarField& p,
    const Foam::volScalarField& T
) const
{
    NotImplemented;
    return uThermo_->ha(p, T);
}


Foam::tmp<Foam::volScalarField::Internal> Foam::ubRhoThermo::ha
(
    const Foam::volScalarField::Internal& p,
    const Foam::volScalarField::Internal& T
) const
{
    NotImplemented;
    return uThermo_->ha(p, T);
}


Foam::tmp<Foam::scalarField> Foam::ubRhoThermo::ha
(
    const Foam::scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return uThermo_->ha(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::ubRhoThermo::ha
(
    const Foam::scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return uThermo_->ha(T, patchi);
}


Foam::tmp<Foam::scalarField> Foam::ubRhoThermo::Cp
(
    const Foam::scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return uThermo_->Cp(T, patchi);
}


Foam::tmp<Foam::scalarField> Foam::ubRhoThermo::Cv
(
    const Foam::scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return uThermo_->Cv(T, patchi);
}


Foam::tmp<Foam::scalarField> Foam::ubRhoThermo::Cpv
(
    const Foam::scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return uThermo_->Cpv(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::ubRhoThermo::The
(
    const Foam::volScalarField& h,
    const Foam::volScalarField& p,
    const Foam::volScalarField& T0
) const
{
    NotImplemented;
    return uThermo_->The(h, p, T0);
}


Foam::tmp<Foam::scalarField> Foam::ubRhoThermo::The
(
    const Foam::scalarField& h,
    const Foam::scalarField& T0,
    const labelList& cells
) const
{
    NotImplemented;
    return uThermo_->The(h, T0, cells);
}


Foam::tmp<Foam::scalarField> Foam::ubRhoThermo::The
(
    const Foam::scalarField& h,
    const Foam::scalarField& T0,
    const label patchi
) const
{
    NotImplemented;
    return uThermo_->The(h, T0, patchi);
}


const Foam::volScalarField& Foam::ubRhoThermo::mu() const
{
    return mu_;
}


const Foam::volScalarField& Foam::ubRhoThermo::kappa() const
{
    return kappa_;
}


// ************************************************************************* //
