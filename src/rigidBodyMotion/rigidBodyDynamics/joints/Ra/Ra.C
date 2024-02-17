/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

#include "Ra.H"
#include "rigidBodyModelState.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace joints
{
    defineTypeNameAndDebug(Ra, 0);

    addToRunTimeSelectionTable
    (
        joint,
        Ra,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::joints::Ra::Ra(const rigidBodyModel& model, const vector& axis)
:
    joint(model, 1)
{
    S_[0] = spatialVector(axis/mag(axis), Zero);
}


Foam::RBD::joints::Ra::Ra(const rigidBodyModel& model, const dictionary& dict)
:
    joint(model, 1)
{
    vector axis(dict.lookup("axis"));
    S_[0] = spatialVector(axis/mag(axis), Zero);
}


Foam::autoPtr<Foam::RBD::joint> Foam::RBD::joints::Ra::clone() const
{
    return autoPtr<joint>(new Ra(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::joints::Ra::~Ra()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::joints::Ra::jcalc
(
    joint::XSvc& J,
    const rigidBodyModelState& state
) const
{
    J.X = Xr(S_[0].w(), state.q()[qIndex_]);
    J.S1 = S_[0];
    J.v = S_[0]*state.qDot()[qIndex_];
    J.c = Zero;
}


void Foam::RBD::joints::Ra::write(Ostream& os) const
{
    joint::write(os);
    writeEntry(os, "axis", S_[0].w());
}


// ************************************************************************* //
