/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Class
    Foam::sixDoFRigidBodyMotionState

Description
    Holds the motion state of sixDoF object.  Wrapped up together
    to allow rapid scatter to other processors.  The processors must all
    maintain exactly the same state data to avoid any drift or inconsistency.

SourceFiles
    sixDoFRigidBodyMotionStateI.H
    sixDoFRigidBodyMotionState.C
    sixDoFRigidBodyMotionStateIO.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyMotionState_H
#define sixDoFRigidBodyMotionState_H

#include "vector.H"
#include "point.H"
#include "diagTensor.H"
#include "tensor.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class sixDoFRigidBodyMotionState;
Istream& operator>>(Istream&, sixDoFRigidBodyMotionState&);
Ostream& operator<<(Ostream&, const sixDoFRigidBodyMotionState&);


/*---------------------------------------------------------------------------*\
                 Class sixDoFRigidBodyMotionState Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyMotionState
{
    // Private Data

        //- Current position of the centre of mass of the body
        point centreOfRotation_;

        //- Orientation, stored as the rotation tensor to transform
        //  from the body to the global reference frame, i.e.:
        //  globalVector = Q_ & bodyLocalVector
        //  bodyLocalVector = Q_.T() & globalVector
        tensor Q_;

        // Linear velocity of body
        vector v_;

        // Total linear acceleration of body
        vector a_;

        //- Angular momentum of body, in body local reference frame
        vector pi_;

        //- Total torque on body, in body local reference frame
        vector tau_;


public:

    // Constructors

        //- Construct null
        sixDoFRigidBodyMotionState();

        //- Construct from dictionary
        sixDoFRigidBodyMotionState(const dictionary& dict);


    //- Destructor
    ~sixDoFRigidBodyMotionState();


    // Member Functions

        // Access

            //- Return access to the centre of mass
            inline const point& centreOfRotation() const;

            //- Return access to the orientation
            inline const tensor& Q() const;

            //- Return access to velocity
            inline const vector& v() const;

            //- Return access to acceleration
            inline const vector& a() const;

            //- Return access to angular momentum
            inline const vector& pi() const;

            //- Return access to torque
            inline const vector& tau() const;


        // Edit

            //- Return non-const access to the centre of mass
            inline point& centreOfRotation();

            //- Return non-const access to the orientation
            inline tensor& Q();

            //- Return non-const access to vector
            inline vector& v();

            //- Return non-const access to acceleration
            inline vector& a();

            //- Return non-const access to angular momentum
            inline vector& pi();

            //- Return non-const access to torque
            inline vector& tau();


        //- Write to dictionary
        void write(dictionary& dict) const;

        //- Write to stream
        void write(Ostream&) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, sixDoFRigidBodyMotionState&);
        friend Ostream& operator<<(Ostream&, const sixDoFRigidBodyMotionState&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sixDoFRigidBodyMotionStateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
