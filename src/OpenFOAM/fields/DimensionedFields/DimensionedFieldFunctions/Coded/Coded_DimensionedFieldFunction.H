/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::DimensionedFieldFunctions::Coded

Description

SourceFiles
    Coded_DimensionedFieldFunction.C

\*---------------------------------------------------------------------------*/

#ifndef Coded_DimensionedFieldFunction_H
#define Coded_DimensionedFieldFunction_H

#include "DimensionedFieldFunction.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace DimensionedFieldFunctions
{

/*---------------------------------------------------------------------------*\
                            Class Coded Declaration
\*---------------------------------------------------------------------------*/

template<class DimensionedFieldType>
class Coded
:
    public DimensionedFieldFunction<DimensionedFieldType>,
    public codedBase
{
    // Private Data

        //- Keywords associated with source code
        static const wordList codeKeys;

        //- Name of the dictionary variables in the source code
        static const wordList codeDictVars;

        // Name of the code options template file
        static const word codeOptions;

        // Names of the source files to compile
        static const wordList compileFiles;

        // Names of the source files to copy (header files)
        static const wordList copyFiles;

        //- The dynamically generated DimensionedFieldFunction pointer
        mutable autoPtr<DimensionedFieldFunction<DimensionedFieldType>>
            redirectFunctionPtr_;


public:

    //- Runtime type information
    TypeName("coded");


    // Constructors

        //- Construct with dictionary to initialise given field
        Coded
        (
            const dictionary& dict,
            DimensionedFieldType& field
        );


    //- Destructor
    virtual ~Coded()
    {}


    // Member Functions

        //- Evaluate the function and set the field
        virtual void evaluate();

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace DimensionedFieldFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "Coded_DimensionedFieldFunction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
