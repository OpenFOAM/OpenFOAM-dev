/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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

#include "solidDisplacementThermo.H"
#include "fvMesh.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(solidDisplacementThermo, 0);
}


void Foam::solidDisplacementThermo::readProperty(volScalarField& prop) const
{
    const dictionary& propDict(subDict(prop.name()));
    const word propType(propDict.lookup("type"));

    if (propType == "uniform")
    {
        prop == dimensionedScalar
        (
            prop.name(),
            prop.dimensions(),
            propDict.lookup<scalar>("value")
        );
    }
    else if (propType == "field")
    {
        const volScalarField propField
        (
            IOobject
            (
                prop.name(),
                prop.mesh().time().timeName(0),
                prop.mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            prop.mesh()
        );
        prop == propField;
    }
    else
    {
        FatalErrorInFunction
            << "Valid type entries are uniform or field for " << prop.name()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidDisplacementThermo::solidDisplacementThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    solidThermo::composite(mesh, phaseName),
    planeStress_(lookup("planeStress")),
    thermalStress_(lookup("thermalStress")),
    Cp_
    (
        IOobject
        (
            phasePropertyName("Cp"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass/dimTemperature
    ),
    kappa_
    (
        IOobject
        (
            phasePropertyName("kappa"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        Cp_.dimensions()*dimensionSet(1, -1, -1, 0, 0)
    ),
    E_
    (
        IOobject
        (
            phasePropertyName("E"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimPressure
    ),
    nu_
    (
        IOobject
        (
            phasePropertyName("nu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    ),
    alphav_
    (
        IOobject
        (
            phasePropertyName("alphav"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimTemperature
    )
{
    readProperty(rho_);
    readProperty(Cp_);
    readProperty(kappa_);
    readProperty(E_);
    readProperty(nu_);
    readProperty(alphav_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidDisplacementThermo::~solidDisplacementThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::solidDisplacementThermo::rho() const
{
    return rho_;
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::rho
(
    const label patchi
) const
{
    return rho_.boundaryField()[patchi];
}


const Foam::volScalarField& Foam::solidDisplacementThermo::E() const
{
    return E_;
}


const Foam::scalarField& Foam::solidDisplacementThermo::E
(
    const label patchi
) const
{
    return E_.boundaryField()[patchi];
}


const Foam::volScalarField& Foam::solidDisplacementThermo::nu() const
{
    return nu_;
}


const Foam::scalarField& Foam::solidDisplacementThermo::nu
(
    const label patchi
) const
{
    return nu_.boundaryField()[patchi];
}


const Foam::volScalarField& Foam::solidDisplacementThermo::alphav() const
{
    return alphav_;
}


const Foam::scalarField& Foam::solidDisplacementThermo::alphav
(
    const label patchi
) const
{
    return alphav_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::solidDisplacementThermo::he()
{
    NotImplemented;
    return rho_;
}


const Foam::volScalarField& Foam::solidDisplacementThermo::he() const
{
    NotImplemented;
    return rho_;
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::he
(
    const scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}

Foam::tmp<Foam::volScalarField> Foam::solidDisplacementThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::he
(
    const scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::solidDisplacementThermo::hs() const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::solidDisplacementThermo::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::hs
(
    const scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::solidDisplacementThermo::ha() const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::solidDisplacementThermo::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::ha
(
    const scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::ha
(
    const scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::solidDisplacementThermo::hc() const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::solidDisplacementThermo::THE
(
    const volScalarField& h,
    const volScalarField& p,
    const volScalarField& T0
) const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::THE
(
    const scalarField& he,
    const scalarField& T0,
    const labelList& cells
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::THE
(
    const scalarField& he,
    const scalarField& T0,
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::solidDisplacementThermo::Cp() const
{
    return Cp_;
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    return Cp_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::solidDisplacementThermo::Cv() const
{
    return Cp_;
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    return Cp_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::solidDisplacementThermo::Cpv() const
{
    return Cp_;
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    return Cp_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::solidDisplacementThermo::CpByCpv() const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::CpByCpv
(
    const scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}



Foam::tmp<Foam::volScalarField> Foam::solidDisplacementThermo::kappa() const
{
    return kappa_;
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::kappa
(
    const label patchi
) const
{
    return kappa_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::solidDisplacementThermo::alphahe() const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::alphahe
(
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::solidDisplacementThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::solidDisplacementThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::volVectorField> Foam::solidDisplacementThermo::Kappa() const
{
    NotImplemented;
    return tmp<volVectorField>(nullptr);
}


Foam::tmp<Foam::vectorField> Foam::solidDisplacementThermo::Kappa
(
    const label patchi
) const
{
    NotImplemented;
    return tmp<vectorField>(nullptr);
}


void Foam::solidDisplacementThermo::correct()
{}


bool Foam::solidDisplacementThermo::read()
{
    return regIOobject::read();
}


// ************************************************************************* //
