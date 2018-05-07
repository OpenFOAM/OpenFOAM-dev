/*---------------------------------------------------------------------------*\3
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

Class
    Foam::RBD::masslessBody

Description

SourceFiles
    masslessBodyI.H
    masslessBody.C

\*---------------------------------------------------------------------------*/

#ifndef RBD_masslessBody_H
#define RBD_masslessBody_H

#include "rigidBody.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{

/*---------------------------------------------------------------------------*\
                         Class masslessBody Declaration
\*---------------------------------------------------------------------------*/

class masslessBody
:
    public rigidBody
{

public:

    //- Runtime type information
    TypeName("masslessBody");


    // Constructors

        //- Construct a massless body
        inline masslessBody();

        //- Construct a named massless body
        inline masslessBody(const word& name);

        //- Construct from dictionary
        inline masslessBody
        (
            const word& name,
            const dictionary& dict
        );

        //- Return clone of this masslessBody
        virtual autoPtr<rigidBody> clone() const;


    //- Destructor
    virtual ~masslessBody();


    // Member Functions

        //- Return true if this body is a massless component of a composite body
        virtual bool massless() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "masslessBodyI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
