/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::functionObjects::div

Group
    grpFieldFunctionObjects

Description
    This function object calculates the divergence of a field.  The operation is
    limited to surfaceScalarFields and volVectorFields, and the output is a
    volScalarField.

SeeAlso
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    div.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_div_H
#define functionObjects_div_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dimensionSet;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class div Declaration
\*---------------------------------------------------------------------------*/

class div
:
    public fvMeshFunctionObject
{
    // Private member data

        //- Name of field to process
        word fieldName_;

        //- Name of result field
        word resultName_;


    // Private Member Functions

        //- Helper function to create/store/return the divergence field
        volScalarField& divField
        (
            const word& gradName,
            const dimensionSet& dims
        );

        //- Helper function to calculate the divergence of different field types
        template<class FieldType>
        void calcDiv
        (
            const word& fieldName,
            const word& resultName,
            bool& processed
        );

        //- Disallow default bitwise copy construct
        div(const div&);

        //- Disallow default bitwise assignment
        void operator=(const div&);


public:

    //- Runtime type information
    TypeName("div");


    // Constructors

        //- Construct from Time and dictionary
        div
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~div();


    // Member Functions

        //- Read the div data
        virtual bool read(const dictionary&);

        //- Calculate the divergence field
        virtual bool execute(const bool postProcess = false);

        //- Write the divergence field
        virtual bool write(const bool postProcess = false);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "divTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
