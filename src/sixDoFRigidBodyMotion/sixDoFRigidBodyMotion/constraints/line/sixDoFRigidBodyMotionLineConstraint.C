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

#include "sixDoFRigidBodyMotionLineConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionConstraints
{
    defineTypeNameAndDebug(line, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionConstraint,
        line,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::line::line
(
    const word& name,
    const dictionary& sDoFRBMCDict,
    const sixDoFRigidBodyMotion& motion
)
:
    sixDoFRigidBodyMotionConstraint(name, sDoFRBMCDict, motion)
{
    read(sDoFRBMCDict);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::line::~line()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionConstraints::line::setCentreOfRotation
(
    point& CofR
) const
{
    CofR = centreOfRotation_;
}


void Foam::sixDoFRigidBodyMotionConstraints::line::constrainTranslation
(
    pointConstraint& pc
) const
{
    pc.combine(pointConstraint(Tuple2<label, vector>(2, direction_)));
}


void Foam::sixDoFRigidBodyMotionConstraints::line::constrainRotation
(
    pointConstraint& pc
) const
{}


bool Foam::sixDoFRigidBodyMotionConstraints::line::read
(
    const dictionary& sDoFRBMCDict
)
{
    sixDoFRigidBodyMotionConstraint::read(sDoFRBMCDict);

    centreOfRotation_ = sDoFRBMCCoeffs_.lookupOrDefault
    (
        "centreOfRotation",
        motion_.initialCentreOfMass()
    );

    sDoFRBMCCoeffs_.lookup("direction") >> direction_;

    scalar magDir(mag(direction_));

    if (magDir > VSMALL)
    {
        direction_ /= magDir;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::sixDoFRigidBodyMotionConstraints::line::read"
            "("
                "const dictionary& sDoFRBMCDict"
            ")"
        )
            << "line direction has zero length"
            << abort(FatalError);
    }

    return true;
}


void Foam::sixDoFRigidBodyMotionConstraints::line::write
(
    Ostream& os
) const
{
    os.writeKeyword("centreOfRotation")
        << centreOfRotation_ << token::END_STATEMENT << nl;
    os.writeKeyword("direction")
        << direction_ << token::END_STATEMENT << nl;
}

// ************************************************************************* //
