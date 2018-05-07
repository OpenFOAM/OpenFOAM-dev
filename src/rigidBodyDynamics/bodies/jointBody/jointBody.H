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
    Foam::RBD::jointBody

Description

SourceFiles
    jointBodyI.H
    jointBody.C

\*---------------------------------------------------------------------------*/

#ifndef RBD_jointBody_H
#define RBD_jointBody_H

#include "masslessBody.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{

/*---------------------------------------------------------------------------*\
                         Class jointBody Declaration
\*---------------------------------------------------------------------------*/

class jointBody
:
    public masslessBody
{

public:

    //- Runtime type information
    TypeName("jointBody");


    // Constructors

        //- Construct a joint body
        inline jointBody();

        //- Construct a named joint body
        inline jointBody(const word& name);

        //- Construct from dictionary
        inline jointBody
        (
            const word& name,
            const dictionary& dict
        );

        //- Return clone of this jointBody
        virtual autoPtr<rigidBody> clone() const;


    //- Destructor
    virtual ~jointBody();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "jointBodyI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
