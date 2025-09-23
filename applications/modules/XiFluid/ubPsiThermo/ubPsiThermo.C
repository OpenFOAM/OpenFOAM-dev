/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "ubPsiThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ubPsiThermo, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ubPsiThermo::ubPsiThermo(const fvMesh& mesh)
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
    uThermo_(ubPsiMulticomponentThermo::New(mesh, unburntPhaseName_)),
    bThermo_(ubPsiMulticomponentThermo::New(mesh, burntPhaseName_)),
    psi_("psi", 1.0/(b_/uThermo_->psi() + c_/bThermo_->psi())),
    mu_("mu", b_*uThermo_->mu() + c_*bThermo_->mu()),
    kappa_("kappa", b_*uThermo_->kappa() + c_*bThermo_->kappa())
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

Foam::ubPsiThermo::~ubPsiThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::IOdictionary& Foam::ubPsiThermo::properties() const
{
    NotImplemented;
    return uThermo_->properties();
}


Foam::IOdictionary& Foam::ubPsiThermo::properties()
{
    NotImplemented;
    return uThermo_->properties();
}


const Foam::fvMesh& Foam::ubPsiThermo::mesh() const
{
    return uThermo_->mesh();
}


const Foam::word& Foam::ubPsiThermo::phaseName() const
{
    return word::null;
}


Foam::word Foam::ubPsiThermo::thermoName() const
{
    return type();
}


void Foam::ubPsiThermo::correct()
{
    uThermo_->correct();
    bThermo_->correct();

    psi_ = 1.0/(b_/uThermo_->psi() + c_/bThermo_->psi());
    mu_ = b_*uThermo_->mu() + c_*bThermo_->mu();
    kappa_ = b_*uThermo_->kappa() + c_*bThermo_->kappa();
}


void Foam::ubPsiThermo::reset()
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


Foam::tmp<Foam::volScalarField> Foam::ubPsiThermo::W() const
{
    NotImplemented;
    return uThermo_->W();
}


Foam::tmp<Foam::scalarField> Foam::ubPsiThermo::W(const label patchi) const
{
    NotImplemented;
    return uThermo_->W(patchi);
}


const Foam::volScalarField& Foam::ubPsiThermo::p() const
{
    return uThermo_->p();
}


Foam::volScalarField& Foam::ubPsiThermo::p()
{
    return uThermo_->p();
}


const Foam::volScalarField& Foam::ubPsiThermo::psi() const
{
    return psi_;
}


const Foam::volScalarField& Foam::ubPsiThermo::T() const
{
    NotImplemented;
    return uThermo_->T();
}


Foam::volScalarField& Foam::ubPsiThermo::T()
{
    NotImplemented;
    return uThermo_->T();
}

const Foam::volScalarField& Foam::ubPsiThermo::he() const
{
    NotImplemented;
    return uThermo_->he();
}

Foam::volScalarField& Foam::ubPsiThermo::he()
{
    NotImplemented;
    return uThermo_->he();
}

const Foam::volScalarField& Foam::ubPsiThermo::Cp() const
{
    NotImplemented;
    return uThermo_->Cp();
}

const Foam::volScalarField& Foam::ubPsiThermo::Cv() const
{
    NotImplemented;
    return uThermo_->Cv();
}

const Foam::volScalarField& Foam::ubPsiThermo::Cpv() const
{
    NotImplemented;
    return uThermo_->Cpv();
}


Foam::tmp<Foam::volScalarField> Foam::ubPsiThermo::rho() const
{
    return p()*psi_;
}


Foam::tmp<Foam::scalarField> Foam::ubPsiThermo::rho(const label patchi) const
{
    return p().boundaryField()[patchi]*psi_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::ubPsiThermo::he
(
    const Foam::volScalarField& p,
    const Foam::volScalarField& T
) const
{
    NotImplemented;
    return uThermo_->he(p, T);
}


Foam::tmp<Foam::volScalarField::Internal> Foam::ubPsiThermo::he
(
    const Foam::volScalarField::Internal& p,
    const Foam::volScalarField::Internal& T
) const
{
    NotImplemented;
    return uThermo_->he(p, T);
}


Foam::tmp<Foam::scalarField> Foam::ubPsiThermo::he
(
    const Foam::scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return uThermo_->he(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::ubPsiThermo::he
(
    const Foam::scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return uThermo_->he(T, patchi);
}


Foam::tmp<Foam::volScalarField::Internal> Foam::ubPsiThermo::he
(
    const Foam::volScalarField::Internal& T,
    const fvSource& model,
    const Foam::volScalarField::Internal& source
) const
{
    NotImplemented;
    return uThermo_->he(T, model, source);
}


Foam::tmp<Foam::scalarField> Foam::ubPsiThermo::he
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


Foam::tmp<Foam::volScalarField> Foam::ubPsiThermo::hs() const
{
    NotImplemented;
    return uThermo_->hs();
}

Foam::tmp<Foam::volScalarField> Foam::ubPsiThermo::hs
(
    const Foam::volScalarField& p,
    const Foam::volScalarField& T
) const
{
    NotImplemented;
    return uThermo_->hs(p, T);
}


Foam::tmp<Foam::volScalarField::Internal> Foam::ubPsiThermo::hs
(
    const Foam::volScalarField::Internal& p,
    const Foam::volScalarField::Internal& T
) const
{
    NotImplemented;
    return uThermo_->hs(p, T);
}


Foam::tmp<Foam::scalarField> Foam::ubPsiThermo::hs
(
    const Foam::scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return uThermo_->hs(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::ubPsiThermo::hs
(
    const Foam::scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return uThermo_->hs(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::ubPsiThermo::ha() const
{
    NotImplemented;
    return uThermo_->ha();
}


Foam::tmp<Foam::volScalarField> Foam::ubPsiThermo::ha
(
    const Foam::volScalarField& p,
    const Foam::volScalarField& T
) const
{
    NotImplemented;
    return uThermo_->ha(p, T);
}


Foam::tmp<Foam::volScalarField::Internal> Foam::ubPsiThermo::ha
(
    const Foam::volScalarField::Internal& p,
    const Foam::volScalarField::Internal& T
) const
{
    NotImplemented;
    return uThermo_->ha(p, T);
}


Foam::tmp<Foam::scalarField> Foam::ubPsiThermo::ha
(
    const Foam::scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return uThermo_->ha(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::ubPsiThermo::ha
(
    const Foam::scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return uThermo_->ha(T, patchi);
}


Foam::tmp<Foam::scalarField> Foam::ubPsiThermo::Cp
(
    const Foam::scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return uThermo_->Cp(T, patchi);
}


Foam::tmp<Foam::scalarField> Foam::ubPsiThermo::Cv
(
    const Foam::scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return uThermo_->Cv(T, patchi);
}


Foam::tmp<Foam::scalarField> Foam::ubPsiThermo::Cpv
(
    const Foam::scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return uThermo_->Cpv(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::ubPsiThermo::The
(
    const Foam::volScalarField& h,
    const Foam::volScalarField& p,
    const Foam::volScalarField& T0
) const
{
    NotImplemented;
    return uThermo_->The(h, p, T0);
}


Foam::tmp<Foam::scalarField> Foam::ubPsiThermo::The
(
    const Foam::scalarField& h,
    const Foam::scalarField& T0,
    const labelList& cells
) const
{
    NotImplemented;
    return uThermo_->The(h, T0, cells);
}


Foam::tmp<Foam::scalarField> Foam::ubPsiThermo::The
(
    const Foam::scalarField& h,
    const Foam::scalarField& T0,
    const label patchi
) const
{
    NotImplemented;
    return uThermo_->The(h, T0, patchi);
}


const Foam::volScalarField& Foam::ubPsiThermo::mu() const
{
    return mu_;
}


const Foam::volScalarField& Foam::ubPsiThermo::kappa() const
{
    return kappa_;
}


// ************************************************************************* //
