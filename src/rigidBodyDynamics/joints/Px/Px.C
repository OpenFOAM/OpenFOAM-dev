/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "Px.H"
#include "rigidBodyModelState.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace joints
{
    defineTypeNameAndDebug(Px, 0);

    addToRunTimeSelectionTable
    (
        joint,
        Px,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::joints::Px::Px(const rigidBodyModel& model)
:
    joint(model, 1)
{
    S_[0] = spatialVector(0, 0, 0, 1, 0, 0);
}


Foam::RBD::joints::Px::Px(const rigidBodyModel& model, const dictionary& dict)
:
    joint(model, 1)
{
    S_[0] = spatialVector(0, 0, 0, 1, 0, 0);
}


Foam::autoPtr<Foam::RBD::joint> Foam::RBD::joints::Px::clone() const
{
    return autoPtr<joint>(new Px(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::joints::Px::~Px()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::joints::Px::jcalc
(
    joint::XSvc& J,
    const rigidBodyModelState& state
) const
{
    J.X = Xt(S_[0].l()*state.q()[qIndex_]);
    J.S1 = S_[0];
    J.v = S_[0]*state.qDot()[qIndex_];
    J.c = Zero;
}


// ************************************************************************* //
