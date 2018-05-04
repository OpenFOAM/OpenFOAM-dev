/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "Pz.H"
#include "rigidBodyModelState.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace joints
{
    defineTypeNameAndDebug(Pz, 0);

    addToRunTimeSelectionTable
    (
        joint,
        Pz,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::joints::Pz::Pz(const rigidBodyModel& model)
:
    joint(model, 1)
{
    S_[0] = spatialVector(0, 0, 0, 0, 0, 1);
}


Foam::RBD::joints::Pz::Pz(const rigidBodyModel& model, const dictionary& dict)
:
    joint(model, 1)
{
    S_[0] = spatialVector(0, 0, 0, 0, 0, 1);
}


Foam::autoPtr<Foam::RBD::joint> Foam::RBD::joints::Pz::clone() const
{
    return autoPtr<joint>(new Pz(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::joints::Pz::~Pz()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::joints::Pz::jcalc
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
