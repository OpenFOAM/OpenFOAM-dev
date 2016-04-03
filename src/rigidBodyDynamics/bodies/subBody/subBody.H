/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::subBody

Description
    This specialized rigidBody holds the original body after it has been merged
    into a parent.

SourceFiles
    subBodyI.H
    subBody.C

\*---------------------------------------------------------------------------*/

#ifndef RBD_subBody_H
#define RBD_subBody_H

#include "rigidBody.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of Foam classes
class Istream;
class Ostream;

namespace RBD
{

// Forward declaration of friend functions and operators
class subBody;
Ostream& operator<<(Ostream&, const subBody&);


/*---------------------------------------------------------------------------*\
                         Class subBody Declaration
\*---------------------------------------------------------------------------*/

class subBody
:
    public rigidBody
{
    // Private data

        //- Original body from which this sub-body was constructed
        autoPtr<rigidBody> body_;

        //- Parent body ID
        label parentID_;

        //- Transform with respect to parent body
        spatialTransform parentXT_;


public:

    // Constructors

        //- Construct a merged version of the given rigidBody
        //  providing the ID of the parent body to which this will be merged
        //  and the transform relative to the parent
        inline subBody
        (
            autoPtr<rigidBody> bodyPtr,
            const label parentID,
            const spatialTransform& parentXT
        );


    // Member Functions

        //- Return the original body from which this sub-body was constructed
        inline const rigidBody& body() const;

        //- Return the parent body Id
        inline label parentID() const;

        //- Return the transform with respect to the parent body
        inline const spatialTransform& parentXT() const;


    // IOstream Operators

        //friend Ostream& operator<<(Ostream&, const subBody&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "subBodyI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
