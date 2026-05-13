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
    Foam::DimensionedFieldFunctions::FieldFunction

Description
    Field function which looks up another field and uses it as the argument to
    evaluate a given Foam::Function1. The user must specify both the name of
    the looked-up field and the function1.

Usage
    To create a temperature dependent electrical conductivity in the
    \c solidElectricalConduction fvModel:
    \verbatim
    solidElectricalConduction
    {
        type        solidElectricalConduction;

        sigma
        {
            type        fieldFunction;
            field       T;
            function
            {
                type        polynomial;
                coeffs      (7.61e7 -1.38e5);
            }
        }
    }
    \endverbatim

SourceFiles
    FieldFunction_DimensionedFieldFunction.C

\*---------------------------------------------------------------------------*/

#ifndef FieldFunction_DimensionedFieldFunction_H
#define FieldFunction_DimensionedFieldFunction_H

#include "DimensionedFieldFunction.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace DimensionedFieldFunctions
{

/*---------------------------------------------------------------------------*\
                            Class FieldFunction Declaration
\*---------------------------------------------------------------------------*/

template<class DimensionedFieldType>
class FieldFunction
:
    public DimensionedFieldFunction<DimensionedFieldType>
{
    // Private Data

        //- Field name
        const word fieldName_;

        //- Field Function1 pointer
        autoPtr<Foam::Function1<typename DimensionedFieldType::Type_>> funcPtr_;

        //- Cached dictionary
        autoPtr<dictionary> dictPtr_;


public:

    //- Runtime type information
    TypeName("fieldFunction");


    // Constructors

        //- Construct with dictionary to initialise given field
        FieldFunction
        (
            const dictionary& dict,
            DimensionedFieldType& field
        );

        //- Construct a copy for the given field
        FieldFunction
        (
            const FieldFunction& dff,
            DimensionedFieldType& field
        );

        //- Construct and return a clone for the specified field
        virtual autoPtr<DimensionedFieldFunction<DimensionedFieldType>>
        clone(DimensionedFieldType& field) const;


    //- Destructor
    virtual ~FieldFunction()
    {}


    // Member Functions

        //- Evaluate the function and set the field
        virtual void evaluate();

        //- Update the field for the current time
        virtual bool update();

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace DimensionedFieldFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "FieldFunction_DimensionedFieldFunction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
