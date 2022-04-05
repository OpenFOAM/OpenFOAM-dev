/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
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

#include "rotating.H"
#include "rigidBodyModelState.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace joints
{
    defineTypeNameAndDebug(rotating, 0);

    addToRunTimeSelectionTable
    (
        joint,
        rotating,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::joints::rotating::rotating
(
    const rigidBodyModel& model,
    const dictionary& dict
)
:
    joint(model, 0),
    omega_(Function1<vector>::New("omega", dict))
{}


Foam::autoPtr<Foam::RBD::joint> Foam::RBD::joints::rotating::clone() const
{
    return autoPtr<joint>(new rotating(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::joints::rotating::~rotating()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::joints::rotating::jcalc
(
    joint::XSvc& J,
    const rigidBodyModelState& state
) const
{
    const scalar t = state.t(), deltaT = state.deltaT();

    if (t < vSmall || deltaT < vSmall)
    {
        return;
    }

    const vector omega = omega_->value(t);
    const vector omegaDot = (omega - omega_->value(t - deltaT))/deltaT;
    const vector theta = omega_->integral(0, t);
    const scalar magTheta = mag(theta);
    const tensor R =
        magTheta > vSmall ? quaternion(theta, magTheta).R() : tensor::I;

    J.X = spatialTransform(R, Zero);
    J.S = Zero;
    J.S1 = Zero;
    J.v = spatialVector(omega, Zero);
    J.c = - spatialVector(omegaDot, Zero);
}


// ************************************************************************* //
