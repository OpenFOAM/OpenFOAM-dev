/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "constSolidThermo.H"
#include "addToRunTimeSelectionTable.H"

/* * * * * * * * * * * * * * * Private Static Data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(constSolidThermo, 0);
    addToRunTimeSelectionTable(basicThermo, constSolidThermo, fvMesh);
    addToRunTimeSelectionTable(solidThermo, constSolidThermo, fvMesh);
}


/* * * * * * * * * * * * * * Protected Member Functions   * * * * * * * * * * */

Foam::volScalarField Foam::constSolidThermo::readProperty
(
    const word& name,
    const dimensionSet& dimensions
) const
{
    const dictionary& propDict(subDict(name));
    const word propType(propDict.lookup("type"));

    if (propType == "uniform")
    {
        return volScalarField
        (
            IOobject
            (
                phasePropertyName(name),
                mesh().time().constant(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar(dimensions, propDict.lookup<scalar>("value"))
        );
    }
    else if (propType == "file")
    {
        return volScalarField
        (
            IOobject
            (
                phasePropertyName(name),
                mesh().time().constant(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh()
        );
    }
    else
    {
        FatalErrorInFunction
            << "Valid type entries are 'uniform' or 'file' for " << name
            << abort(FatalError);

        return volScalarField::null();
    }
}


void Foam::constSolidThermo::readProperty
(
    const word& name,
    volScalarField& prop
) const
{
    const dictionary& propDict(subDict(name));
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
    else if (propType == "file")
    {
        const volScalarField propField
        (
            IOobject
            (
                prop.name(),
                prop.mesh().time().constant(),
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
            << "Valid type entries are 'uniform' or 'file' for " << prop.name()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constSolidThermo::constSolidThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    solidThermo::composite(mesh, phaseName),
    Cv_(readProperty("Cv", dimEnergy/dimMass/dimTemperature)),
    e_
    (
        IOobject
        (
            phasePropertyName("e"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Cv_*T_,
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    )
{
    readProperty("rho", rho_);
    readProperty("kappa", kappa_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constSolidThermo::~constSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::constSolidThermo::Cpv() const
{
    return Cv_;
}


Foam::volScalarField& Foam::constSolidThermo::he()
{
    return e_;
}


const Foam::volScalarField& Foam::constSolidThermo::he() const
{
    return e_;
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::he
(
    const scalarField& T,
    const labelList& cells
) const
{
    return scalarField(Cv_, cells)*T;
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return Cv_*T;
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::he
(
    const scalarField& T,
    const label patchi
) const
{
    return Cv_.boundaryField()[patchi]*T;
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermo::hs() const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermo::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::hs
(
    const scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermo::ha() const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermo::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::ha
(
    const scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::ha
(
    const scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermo::hc() const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermo::THE
(
    const volScalarField& h,
    const volScalarField& p,
    const volScalarField& T0
) const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::THE
(
    const scalarField& he,
    const scalarField& T0,
    const labelList& cells
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::THE
(
    const scalarField& he,
    const scalarField& T0,
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


const Foam::volScalarField& Foam::constSolidThermo::Cp() const
{
    return Cv_;
}


const Foam::volScalarField& Foam::constSolidThermo::Cv() const
{
    return Cv_;
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    return Cv_.boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    return Cv_.boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    return Cv_.boundaryField()[patchi];
}


Foam::tmp<Foam::volVectorField> Foam::constSolidThermo::Kappa() const
{
    NotImplemented;
    return tmp<volVectorField>(nullptr);
}


Foam::tmp<Foam::vectorField> Foam::constSolidThermo::Kappa
(
    const label patchi
) const
{
    NotImplemented;
    return tmp<vectorField>(nullptr);
}


void Foam::constSolidThermo::correct()
{
    T_ = e_/Cv_;
}


bool Foam::constSolidThermo::read()
{
    return regIOobject::read();
}


// ************************************************************************* //
