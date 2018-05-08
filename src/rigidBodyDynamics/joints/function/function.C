/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "function.H"
#include "rigidBodyModelState.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace joints
{
    defineTypeNameAndDebug(function, 0);

    addToRunTimeSelectionTable
    (
        joint,
        function,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::joints::function::function
(
    const rigidBodyModel& model,
    const dictionary& dict
)
:
    joint(model, 0),
    f_(Function1<scalar>::New("function", dict))
{}


Foam::autoPtr<Foam::RBD::joint> Foam::RBD::joints::function::clone() const
{
    return autoPtr<joint>(new function(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::joints::function::~function()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::joints::function::jcalc
(
    joint::XSvc& J,
    const rigidBodyModelState& state
) const
{
    const label lambda = model_.lambda()[index_];
    const joint& parent = model_.joints()[lambda];

    spatialVector x(Zero), v(Zero), c(Zero);
    for(label i = 0; i < model_.joints()[lambda].nDoF(); ++ i)
    {
        x += f_->value(state.q()[parent.qIndex() + i])*parent.S()[i];
        v += f_->value(state.qDot()[parent.qIndex() + i])*parent.S()[i];
        c += f_->value(state.qDdot()[parent.qIndex() + i])*parent.S()[i];
    }

    const scalar magW = mag(x.w());

    const tensor X(magW > vSmall ? quaternion(x.w(), magW).R() : tensor::I);

    J.X = spatialTransform(X, x.l());
    J.S = Zero;
    J.S1 = Zero;
    J.v = v;
    J.c = c;
}


// ************************************************************************* //
