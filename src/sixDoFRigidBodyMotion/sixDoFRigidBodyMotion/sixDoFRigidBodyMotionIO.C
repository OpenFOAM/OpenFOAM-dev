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

#include "sixDoFRigidBodyMotion.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::sixDoFRigidBodyMotion::read(const dictionary& dict)
{
    dict.lookup("mass") >> mass_;
    dict.lookup("momentOfInertia") >> momentOfInertia_;
    aRelax_ = dict.lookupOrDefault<scalar>("accelerationRelaxation", 1.0);
    aDamp_ = dict.lookupOrDefault<scalar>("accelerationDamping", 1.0);
    report_ = dict.lookupOrDefault<Switch>("report", false);

    restraints_.clear();
    addRestraints(dict);

    constraints_.clear();
    addConstraints(dict);

    return true;
}


void Foam::sixDoFRigidBodyMotion::write(Ostream& os) const
{
    motionState_.write(os);

    os.writeKeyword("centreOfMass")
        << initialCentreOfMass_ << token::END_STATEMENT << nl;
    os.writeKeyword("initialOrientation")
        << initialQ_ << token::END_STATEMENT << nl;
    os.writeKeyword("mass")
        << mass_ << token::END_STATEMENT << nl;
    os.writeKeyword("momentOfInertia")
        << momentOfInertia_ << token::END_STATEMENT << nl;
    os.writeKeyword("accelerationRelaxation")
        << aRelax_ << token::END_STATEMENT << nl;
    os.writeKeyword("accelerationDamping")
        << aDamp_ << token::END_STATEMENT << nl;
    os.writeKeyword("report")
        << report_ << token::END_STATEMENT << nl;

    if (!restraints_.empty())
    {
        os  << indent << "restraints" << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;

        forAll(restraints_, rI)
        {
            word restraintType = restraints_[rI].type();

            os  << indent << restraints_[rI].name() << nl
                << indent << token::BEGIN_BLOCK << incrIndent << endl;

            os.writeKeyword("sixDoFRigidBodyMotionRestraint")
                << restraintType << token::END_STATEMENT << nl;

            restraints_[rI].write(os);

            os  << decrIndent << indent << token::END_BLOCK << endl;
        }

        os  << decrIndent << indent << token::END_BLOCK << nl;
    }

    if (!constraints_.empty())
    {
        os  << indent << "constraints" << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;

        forAll(constraints_, rI)
        {
            word constraintType = constraints_[rI].type();

            os  << indent << constraints_[rI].name() << nl
                << indent << token::BEGIN_BLOCK << incrIndent << endl;

            os.writeKeyword("sixDoFRigidBodyMotionConstraint")
                << constraintType << token::END_STATEMENT << nl;

            constraints_[rI].sixDoFRigidBodyMotionConstraint::write(os);

            constraints_[rI].write(os);

            os  << decrIndent << indent << token::END_BLOCK << endl;
        }

        os  << decrIndent << indent << token::END_BLOCK << nl;
    }
}


// ************************************************************************* //
