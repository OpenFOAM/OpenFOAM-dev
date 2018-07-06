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

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

inline Foam::RBD::rigidBodyModelState& Foam::RBD::rigidBodySolver::state()
{
    return model_.motionState_;
}


inline Foam::scalarField& Foam::RBD::rigidBodySolver::q()
{
    return state().q();
}


inline Foam::scalarField& Foam::RBD::rigidBodySolver::qDot()
{
    return state().qDot();
}


inline Foam::scalarField& Foam::RBD::rigidBodySolver::qDdot()
{
    return state().qDdot();
}


inline Foam::scalar Foam::RBD::rigidBodySolver::deltaT() const
{
    return model_.motionState_.deltaT();
}


inline const Foam::RBD::rigidBodyModelState&
Foam::RBD::rigidBodySolver::state0() const
{
    return model_.motionState0_;
}

inline const Foam::scalarField& Foam::RBD::rigidBodySolver::q0() const
{
    return state0().q();
}


inline const Foam::scalarField& Foam::RBD::rigidBodySolver::qDot0() const
{
    return state0().qDot();
}


inline const Foam::scalarField& Foam::RBD::rigidBodySolver::qDdot0() const
{
    return state0().qDdot();
}


inline Foam::scalar Foam::RBD::rigidBodySolver::deltaT0() const
{
    return model_.motionState0_.deltaT();
}


// ************************************************************************* //
