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

#include "spherical.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(spherical, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::clouds::spherical::dIo(const cloud& c)
{
    return
        IOobject
        (
            "d",
            c.time().name(),
            c.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::clouds::spherical::calcv
(
    const LagrangianModelRef& model,
    const LagrangianSubMesh& subMesh
) const
{
    return constant::mathematical::pi*pow3(d(model, subMesh))/6;
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::clouds::spherical::calca
(
    const LagrangianModelRef& model,
    const LagrangianSubMesh& subMesh
) const
{
    return constant::mathematical::pi*sqr(d(model, subMesh));
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::clouds::spherical::correct(const LagrangianSubScalarSubField& v)
{
    const LagrangianSubMesh& subMesh = v.mesh();

    LagrangianSubScalarSubField d(subMesh.sub(this->d));

    d = cbrt(6/constant::mathematical::pi*v);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::spherical::spherical(const cloud& c)
:
    shaped(c),
    d(c.stateField<scalar>(dIo(c), c.mesh()))
{}


Foam::clouds::spherical::spherical(const cloud& c, const grouped& groupedCloud)
:
    shaped(c, groupedCloud),
    d(c.stateField<scalar>(dIo(c), c.mesh()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::spherical::~spherical()
{}


// ************************************************************************* //
