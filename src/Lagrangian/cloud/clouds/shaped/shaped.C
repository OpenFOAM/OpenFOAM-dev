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

#include "shaped.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(shaped, 0);
}
}


const Foam::word Foam::clouds::shaped::vName("v");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::shaped::shaped(const cloud& c)
:
    v(c.derivedField<scalar>(vName, *this, &shaped::calcv)),
    a(c.derivedField<scalar>(*this, &shaped::calca)),
    alpha
    (
        c.averageField<scalar>
        (
            "alpha",
            refCast<const fvMesh>(c.mesh().mesh()).V(),
            v
        )
    )
{}


Foam::clouds::shaped::shaped(const cloud& c, const grouped& groupedCloud)
:
    v(c.derivedField<scalar>(vName, *this, &shaped::calcv)),
    a(c.derivedField<scalar>(*this, &shaped::calca)),
    alpha
    (
        c.averageField<scalar>
        (
            "alpha",
            refCast<const fvMesh>(c.mesh().mesh()).V(),
            c.derivedField<scalar>
            (
                [&]
                (
                    const LagrangianModelRef& model,
                    const LagrangianSubMesh& subMesh
                )
                {
                    return
                        groupedCloud.number(model, subMesh)
                       *v(model, subMesh);
                }
            )
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::shaped::~shaped()
{}


// ************************************************************************* //
