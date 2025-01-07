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

#include "massive.H"
#include "calculatedLagrangianPatchFields.H"
#include "massLagrangianScalarFieldSource.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(massive, 0);
}
}


const Foam::word Foam::clouds::massive::mName("m");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::massive::massive(const cloud& c, const shaped& shapedCloud)
:
    rho
    (
        c.stateField<scalar>
        (
            IOobject
            (
                "rho",
                c.time().name(),
                c.mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            c.mesh()
        )
    ),
    m
    (
        c.derivedField<scalar>
        (
            mName,
            [&]
            (
                const LagrangianModelRef& model,
                const LagrangianSubMesh& subMesh
            )
            {
                return shapedCloud.v(model, subMesh)*rho(model, subMesh);
            }
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::massive::~massive()
{}


// ************************************************************************* //
