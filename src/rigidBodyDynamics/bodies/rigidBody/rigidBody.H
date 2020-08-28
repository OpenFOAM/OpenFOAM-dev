/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2020 OpenFOAM Foundation
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
    Foam::RBD::rigidBody

Description

SourceFiles
    rigidBodyI.H

\*---------------------------------------------------------------------------*/

#ifndef RBD_rigidBody_H
#define RBD_rigidBody_H

#include "rigidBodyInertia.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{

// Forward declaration of classes
class subBody;

// Forward declaration of friend functions and operators
class rigidBody;

inline Ostream& operator<<(Ostream&, const rigidBody&);


/*---------------------------------------------------------------------------*\
                          Class rigidBody Declaration
\*---------------------------------------------------------------------------*/

class rigidBody
:
    public rigidBodyInertia
{
    // Private Data

        //- Name of body
        word name_;


public:

    //- Runtime type information
    TypeName("rigidBody");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            rigidBody,
            dictionary,
            (const word& name, const dictionary& dict),
            (name, dict)
        );


    // Constructors

        //- Construct from mass, centre of mass and moment of inertia tensor
        //  about the centre of mass
        inline rigidBody
        (
            const word& name,
            const scalar& m,
            const vector& c,
            const symmTensor& Ic
        );

        //- Construct from with the given rigidBodyInertia
        inline rigidBody
        (
            const word& name,
            const rigidBodyInertia& rbi
        );

        //- Construct from dictionary
        inline rigidBody
        (
            const word& name,
            const dictionary& dict
        );

        //- Return clone of this rigidBody
        virtual autoPtr<rigidBody> clone() const;


    // Selectors

        //- Select constructed from components
        static autoPtr<rigidBody> New
        (
            const word& name,
            const scalar& m,
            const vector& c,
            const symmTensor& Ic
        );

        //- Select constructed from name and dictionary
        static autoPtr<rigidBody> New
        (
            const word& name,
            const dictionary& dict
        );


    //- Destructor
    virtual ~rigidBody();


    // Member Functions

        //- Return name
        inline const word& name() const;

        //- Return false as this body is not massless
        virtual bool massless() const;

        //- Merge a body into this parent body
        void merge(const subBody&);

        //- Write
        virtual void write(Ostream&) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const rigidBody&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rigidBodyI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
