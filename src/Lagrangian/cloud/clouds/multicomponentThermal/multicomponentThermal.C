/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "multicomponentThermal.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(multicomponentThermal, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::multicomponentThermal::multicomponentThermal
(
    const cloud& c,
    const shaped& shapedCloud,
    const carried& carriedCloud
)
:
    Thermal<multicomponentLagrangianThermo>(c, shapedCloud, carriedCloud),
    Y()
{
    multicomponentLagrangianThermo& thermo =
        this->thermo<multicomponentLagrangianThermo>();

    Y.resize(thermo.Y().size());

    forAll(Y, i)
    {
        Y.set
        (
            i,
            new CloudStateFieldRef<scalar>
            (
                thermo.Y()[i]
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::multicomponentThermal::~multicomponentThermal()
{}


// ************************************************************************* //
