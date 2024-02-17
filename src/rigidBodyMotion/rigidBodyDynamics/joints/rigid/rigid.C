/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2024 OpenFOAM Foundation
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

#include "rigid.H"
#include "rigidBodyModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace joints
{
    defineTypeNameAndDebug(rigid, 0);

    addToRunTimeSelectionTable
    (
        joint,
        rigid,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::joints::rigid::rigid(const rigidBodyModel& model)
:
    joint(model, 0)
{}


Foam::RBD::joints::rigid::rigid
(
    const rigidBodyModel& model,
    const dictionary& dict
)
:
    joint(model, 0)
{}


Foam::autoPtr<Foam::RBD::joint> Foam::RBD::joints::rigid::clone() const
{
    return autoPtr<joint>(new rigid(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::joints::rigid::~rigid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::joints::rigid::jcalc
(
    joint::XSvc& J,
    const rigidBodyModelState& state
) const
{
    J.X = spatialTransform();
    J.S = Zero;
    J.S1 = Zero;
    J.v = Zero;
    J.c = Zero;
}


// ************************************************************************* //
