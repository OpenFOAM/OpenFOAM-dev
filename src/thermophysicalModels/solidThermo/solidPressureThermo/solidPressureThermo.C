/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "solidPressureThermo.H"
#include "fvMesh.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(solidPressureThermo, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidPressureThermo::solidPressureThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    solidThermo(mesh, phaseName),
    p_(lookupOrConstruct(mesh, "p"))
{}


Foam::solidPressureThermo::solidPressureThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    solidThermo(mesh, dict, phaseName),
    p_(lookupOrConstruct(mesh, "p"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidPressureThermo::~solidPressureThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::solidPressureThermo::p() const
{
    return p_;
}


Foam::volScalarField& Foam::solidPressureThermo::p()
{
    return p_;
}


// ************************************************************************* //
