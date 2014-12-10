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

#include "sixDoFRigidBodyMotionAxisConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionConstraints
{
    defineTypeNameAndDebug(axis, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionConstraint,
        axis,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::axis::axis
(
    const word& name,
    const dictionary& sDoFRBMCDict,
    const sixDoFRigidBodyMotion& motion
)
:
    sixDoFRigidBodyMotionConstraint(name, sDoFRBMCDict, motion),
    axis_()
{
    read(sDoFRBMCDict);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::axis::~axis()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionConstraints::axis::constrainTranslation
(
    pointConstraint& pc
) const
{}


void Foam::sixDoFRigidBodyMotionConstraints::axis::constrainRotation
(
    pointConstraint& pc
) const
{
    pc.combine(pointConstraint(Tuple2<label, vector>(2, axis_)));
}


bool Foam::sixDoFRigidBodyMotionConstraints::axis::read
(
    const dictionary& sDoFRBMCDict
)
{
    sixDoFRigidBodyMotionConstraint::read(sDoFRBMCDict);

    sDoFRBMCCoeffs_.lookup("axis") >> axis_;

    scalar magFixedAxis(mag(axis_));

    if (magFixedAxis > VSMALL)
    {
        axis_ /= magFixedAxis;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::sixDoFRigidBodyMotionConstraints::axis::read"
            "("
                "const dictionary& sDoFRBMCDict"
            ")"
        )
            << "axis has zero length"
            << abort(FatalError);
    }

    return true;
}


void Foam::sixDoFRigidBodyMotionConstraints::axis::write
(
    Ostream& os
) const
{
    os.writeKeyword("axis")
        << axis_ << token::END_STATEMENT << nl;
}

// ************************************************************************* //
