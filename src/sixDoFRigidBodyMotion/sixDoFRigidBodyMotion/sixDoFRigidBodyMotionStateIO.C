/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "sixDoFRigidBodyMotionState.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionState::write(dictionary& dict) const
{
    dict.add("centreOfRotation", centreOfRotation_);
    dict.add("orientation", Q_);
    dict.add("velocity", v_);
    dict.add("acceleration", a_);
    dict.add("angularMomentum", pi_);
    dict.add("torque", tau_);
}


void Foam::sixDoFRigidBodyMotionState::write(Ostream& os) const
{
    os.writeKeyword("centreOfRotation")
        << centreOfRotation_ << token::END_STATEMENT << nl;
    os.writeKeyword("orientation")
        << Q_ << token::END_STATEMENT << nl;
    os.writeKeyword("velocity")
        << v_ << token::END_STATEMENT << nl;
    os.writeKeyword("acceleration")
        << a_ << token::END_STATEMENT << nl;
    os.writeKeyword("angularMomentum")
        << pi_ << token::END_STATEMENT << nl;
    os.writeKeyword("torque")
        << tau_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is, sixDoFRigidBodyMotionState& sDoFRBMS
)
{
    is  >> sDoFRBMS.centreOfRotation_
        >> sDoFRBMS.Q_
        >> sDoFRBMS.v_
        >> sDoFRBMS.a_
        >> sDoFRBMS.pi_
        >> sDoFRBMS.tau_;

    // Check state of Istream
    is.check
    (
        "Foam::Istream& Foam::operator>>"
        "(Foam::Istream&, Foam::sixDoFRigidBodyMotionState&)"
    );

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const sixDoFRigidBodyMotionState& sDoFRBMS
)
{
    os  << token::SPACE << sDoFRBMS.centreOfRotation()
        << token::SPACE << sDoFRBMS.Q()
        << token::SPACE << sDoFRBMS.v()
        << token::SPACE << sDoFRBMS.a()
        << token::SPACE << sDoFRBMS.pi()
        << token::SPACE << sDoFRBMS.tau();

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<(Foam::Ostream&, "
        "const Foam::sixDoFRigidBodyMotionState&)"
    );

    return os;
}


// ************************************************************************* //
