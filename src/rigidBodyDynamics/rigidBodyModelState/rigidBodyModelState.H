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

Class
    Foam::RBD::rigidBodyModelState

Description
    Holds the motion state of rigid-body model.

SourceFiles
    rigidBodyModelStateI.H
    rigidBodyModelState.C
    rigidBodyModelStateIO.C

\*---------------------------------------------------------------------------*/

#ifndef rigidBodyModelState_H
#define rigidBodyModelState_H

#include "rigidBodyModel.H"
#include "scalarField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

namespace RBD
{

// Forward declaration of friend functions and operators
class rigidBodyModelState;
Istream& operator>>(Istream&, rigidBodyModelState&);
Ostream& operator<<(Ostream&, const rigidBodyModelState&);


/*---------------------------------------------------------------------------*\
                  Class rigidBodyModelState Declaration
\*---------------------------------------------------------------------------*/

class rigidBodyModelState
{
    // Private data

        //- Joint position and orientation
        scalarField q_;

        //- Joint velocity
        scalarField qDot_;

        //- Joint acceleration
        scalarField qDdot_;

        //- The time
        scalar t_;

        //- The time-step used to integrate to this state
        scalar deltaT_;


public:

    // Constructors

        //- Construct for the given rigidBodyModel
        rigidBodyModelState(const rigidBodyModel& model);

        //- Construct from dictionary for the given rigidBodyModel
        rigidBodyModelState
        (
            const rigidBodyModel& model,
            const dictionary& dict
        );


    // Member Functions

        // Access

            //- Return access to the joint position and orientation
            inline const scalarField& q() const;

            //- Return access to the joint velocity
            inline const scalarField& qDot() const;

            //- Return access to the joint acceleration
            inline const scalarField& qDdot() const;

            //- Return access to the time
            inline scalar t() const;

            //- Return access to the time-step
            inline scalar deltaT() const;


        // Edit

            //- Return access to the joint position and orientation
            inline scalarField& q();

            //- Return access to the joint velocity
            inline scalarField& qDot();

            //- Return access to the joint acceleration
            inline scalarField& qDdot();

            //- Return access to the time
            inline scalar& t();

            //- Return access to the time-step
            inline scalar& deltaT();


        //- Write to dictionary
        void write(dictionary& dict) const;

        //- Write to stream
        void write(Ostream&) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, rigidBodyModelState&);
        friend Ostream& operator<<(Ostream&, const rigidBodyModelState&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rigidBodyModelStateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
