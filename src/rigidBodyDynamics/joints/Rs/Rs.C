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

#include "Rs.H"
#include "rigidBodyModelState.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace joints
{
    defineTypeNameAndDebug(Rs, 0);

    addToRunTimeSelectionTable
    (
        joint,
        Rs,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::joints::Rs::Rs(const rigidBodyModel& model)
:
    joint(model, 3)
{
    S_[0] = spatialVector(1, 0, 0, 0, 0, 0);
    S_[1] = spatialVector (0, 1, 0, 0, 0, 0);
    S_[2] = spatialVector(0, 0, 1, 0, 0, 0);
}


Foam::RBD::joints::Rs::Rs(const rigidBodyModel& model, const dictionary& dict)
:
    joint(model, 3)
{
    S_[0] = spatialVector(1, 0, 0, 0, 0, 0);
    S_[1] = spatialVector (0, 1, 0, 0, 0, 0);
    S_[2] = spatialVector(0, 0, 1, 0, 0, 0);
}


Foam::autoPtr<Foam::RBD::joint> Foam::RBD::joints::Rs::clone() const
{
    return autoPtr<joint>(new Rs(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::joints::Rs::~Rs()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::RBD::joints::Rs::unitQuaternion() const
{
    return true;
}


void Foam::RBD::joints::Rs::jcalc
(
    joint::XSvc& J,
    const rigidBodyModelState& state
) const
{
    J.X.E() = joint::unitQuaternion(state.q()).R().T();
    J.X.r() = Zero;

    J.S = Zero;
    J.S.xx() = 1;
    J.S.yy() = 1;
    J.S.zz() = 1;

    J.v = spatialVector(state.qDot().block<vector>(qIndex_), Zero);
    J.c = Zero;
}


// ************************************************************************* //
