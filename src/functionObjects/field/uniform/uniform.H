/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::functionObjects::uniform

Description
    Generate a uniform field

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    uniform.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_uniform_H
#define functionObjects_uniform_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class uniform Declaration
\*---------------------------------------------------------------------------*/

class uniform
:
    public fvMeshFunctionObject
{
    // Private Data

        //- The type of the field to be written
        word fieldType_;

        //- The field's name
        word name_;

        //- The field's dimensions
        dimensionSet dimensions_;

        //- The field's value
        #define valueType(Type, none) Type Type##Value_;
        FOR_ALL_FIELD_TYPES(valueType);
        #undef valueType


public:

    //- Runtime type information
    TypeName("uniform");


    // Constructors

        //- Construct from Time and dictionary
        uniform
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    //- Destructor
    virtual ~uniform();


    // Member Functions

        //- Read the uniform data
        virtual bool read(const dictionary&);

        //- Calculate the field
        virtual bool execute();

        //- Write the field
        virtual bool write();

        //- Clear the field from the objectRegistry
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
