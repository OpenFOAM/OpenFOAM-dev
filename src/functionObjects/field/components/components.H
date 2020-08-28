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
    Foam::functionObjects::components

Description
    Calculates the components of a field.

    The operation can be applied to any volume or surface fields generating a
    volume or surface scalar fields for each component.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    components.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_components_H
#define functionObjects_components_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class components Declaration
\*---------------------------------------------------------------------------*/

class components
:
    public fieldExpression
{
    // Private member data

        //- List of the component field names
        wordList resultNames_;


    // Private Member Functions

        //- Calculate the components of the field with the specified type
        //  and register the result
        template<class GeoFieldType>
        bool calcFieldComponents();

        //- Calculate the components of the field with the specified
        //  element type and register the result
        template<class Type>
        bool calcComponents();

        //- Calculate the components of the field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("components");


    // Constructors

        //- Construct from Time and dictionary
        components
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~components();


    // Member Functions

        //- Write the component fields
        virtual bool write();

        //- Clear the component fields from the objectRegistry
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "componentsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
