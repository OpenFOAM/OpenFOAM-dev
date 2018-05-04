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

#include "Rxyz.H"
#include "rigidBodyModelState.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace joints
{
    defineTypeNameAndDebug(Rxyz, 0);

    addToRunTimeSelectionTable
    (
        joint,
        Rxyz,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::joints::Rxyz::Rxyz(const rigidBodyModel& model)
:
    joint(model, 3)
{
    S_[0] = spatialVector(1, 0, 0, 0, 0, 0);
    S_[1] = spatialVector(0, 1, 0, 0, 0, 0);
    S_[2] = spatialVector(0, 0, 1, 0, 0, 0);
}


Foam::RBD::joints::Rxyz::Rxyz
(
    const rigidBodyModel& model,
    const dictionary& dict
)
:
    joint(model, 3)
{
    S_[0] = spatialVector(1, 0, 0, 0, 0, 0);
    S_[1] = spatialVector(0, 1, 0, 0, 0, 0);
    S_[2] = spatialVector(0, 0, 1, 0, 0, 0);
}


Foam::autoPtr<Foam::RBD::joint> Foam::RBD::joints::Rxyz::clone() const
{
    return autoPtr<joint>(new Rxyz(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::joints::Rxyz::~Rxyz()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::joints::Rxyz::jcalc
(
    joint::XSvc& J,
    const rigidBodyModelState& state
) const
{
    vector qj(state.q().block<vector>(qIndex_));

    scalar s0 = sin(qj.x());
    scalar c0 = cos(qj.x());
    scalar s1 = sin(qj.y());
    scalar c1 = cos(qj.y());
    scalar s2 = sin(qj.z());
    scalar c2 = cos(qj.z());

    J.X.E() = tensor
    (
        c2*c1,            s2*c0 + c2*s1*s0, s2*s0 - c2*s1*c0,
       -s2*c1,            c2*c0 - s2*s1*s0, c2*s0 + s2*s1*c0,
        s1,              -c1*s0,            c1*c0
    );
    J.X.r() = Zero;

    J.S = Zero;
    J.S.xx() = c2*c1;
    J.S.xy() = s2;
    J.S.yx() = -s2*c1;
    J.S.yy() = c2;
    J.S.zx() = s1;
    J.S.zz() = 1;

    vector qDotj(state.qDot().block<vector>(qIndex_));
    J.v = J.S & qDotj;

    J.c = spatialVector
    (
       -s2*c1*qDotj.z()*qDotj.x()
      - c2*s1*qDotj.y()*qDotj.x()
      + c2*qDotj.z()*qDotj.y(),

       -c2*c1*qDotj.z()*qDotj.x()
      + s2*s1*qDotj.y()*qDotj.x()
      - s2*qDotj.z()*qDotj.y(),

        c1*qDotj.y()*qDotj.x(),

        0,
        0,
        0
    );
}


// ************************************************************************* //
