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

#include "Pxyz.H"
#include "rigidBodyModelState.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace joints
{
    defineTypeNameAndDebug(Pxyz, 0);

    addToRunTimeSelectionTable
    (
        joint,
        Pxyz,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::joints::Pxyz::Pxyz(const rigidBodyModel& model)
:
    joint(model, 3)
{
    S_[0] = spatialVector(0, 0, 0, 1, 0, 0);
    S_[1] = spatialVector(0, 0, 0, 0, 1, 0);
    S_[2] = spatialVector(0, 0, 0, 0, 0, 1);
}


Foam::RBD::joints::Pxyz::Pxyz
(
    const rigidBodyModel& model,
    const dictionary& dict
)
:
    joint(model, 3)
{
    S_[0] = spatialVector(0, 0, 0, 1, 0, 0);
    S_[1] = spatialVector(0, 0, 0, 0, 1, 0);
    S_[2] = spatialVector(0, 0, 0, 0, 0, 1);
}


Foam::autoPtr<Foam::RBD::joint> Foam::RBD::joints::Pxyz::clone() const
{
    return autoPtr<joint>(new Pxyz(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::joints::Pxyz::~Pxyz()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::joints::Pxyz::jcalc
(
    joint::XSvc& J,
    const rigidBodyModelState& state
) const
{
    J.X.E() = tensor::I;
    J.X.r() = state.q().block<vector>(qIndex_);

    J.S = Zero;
    J.S(3,0) = 1;
    J.S(4,1) = 1;
    J.S(5,2) = 1;

    J.v = spatialVector(Zero, state.qDot().block<vector>(qIndex_));
    J.c = Zero;
}


// ************************************************************************* //
