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
    Foam::DimensionedFieldFunction

Description

SourceFiles
    DimensionedFieldFunction.C

\*---------------------------------------------------------------------------*/

#ifndef DimensionedFieldFunction_H
#define DimensionedFieldFunction_H

#include "DimensionedFieldFwd.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class DimensionedFieldFunction Declaration
\*---------------------------------------------------------------------------*/

template<class DimensionedFieldType>
class DimensionedFieldFunction
{

protected:

    //- Reference to the field the function applies to
    DimensionedFieldType& field_;


public:

    //- Runtime type information
    TypeName("DimensionedFieldFunction");


    // Declare runtime constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            DimensionedFieldFunction,
            dictionary,
            (
                const dictionary& dict,
                DimensionedFieldType& field
            ),
            (dict, field)
        );


    // Constructors

        //- Construct with dictionary to initialise given field
        DimensionedFieldFunction
        (
            const dictionary& dict,
            DimensionedFieldType& field
        );

        //- Construct a copy for the given field
        DimensionedFieldFunction
        (
            const DimensionedFieldFunction& dff,
            DimensionedFieldType& field
        );

        //- Disallow default bitwise copy construction
        DimensionedFieldFunction
        (
            const DimensionedFieldFunction<DimensionedFieldType>&
        ) = delete;

        //- Construct and return a clone for the specified field
        virtual autoPtr<DimensionedFieldFunction<DimensionedFieldType>>
        clone(DimensionedFieldType& field) const = 0;

        //- Construct and return a clone for the specified field
        virtual autoPtr<DimensionedFieldFunction<DimensionedFieldType>>
        clone() const
        {
            NotImplemented;
            return autoPtr<DimensionedFieldFunction<DimensionedFieldType>>
            (
                nullptr
            );
        }


    // Selectors

        //- Select null constructed
        static autoPtr<DimensionedFieldFunction<DimensionedFieldType>> New
        (
            const dictionary& dict,
            DimensionedFieldType& field
        );


    //- Destructor
    virtual ~DimensionedFieldFunction()
    {}


    // Member Functions

        //- Evaluate the function and set the field
        virtual void evaluate() = 0;

        //- Update the field if the function has changed
        virtual void update()
        {}

        //- Reset the field size
        virtual void reset()
        {
            field_.setSize(field_.mesh().size());
            evaluate();
        }

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const DimensionedFieldFunction<DimensionedFieldType>&
        ) = delete;
};


template<class DimensionedFieldType>
void writeEntry
(
    Ostream& os,
    const word& key,
    const DimensionedFieldFunction<DimensionedFieldType>& f
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineDimensionedFieldFunction(DimensionedFieldType)                   \
                                                                               \
    typedef DimensionedFieldFunction<DimensionedFieldType>                     \
        DimensionedFieldType##DimensionedFieldFunction;                        \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        DimensionedFieldType##DimensionedFieldFunction,                        \
        0                                                                      \
    );                                                                         \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        DimensionedFieldType##DimensionedFieldFunction,                        \
        dictionary                                                             \
    );


#define addDimensionedFieldFunction(SS, DimensionedFieldType)                  \
                                                                               \
    typedef DimensionedFieldFunction<DimensionedFieldType>                     \
        DimensionedFieldType##DimensionedFieldFunction;                        \
    typedef SS<DimensionedFieldType>                                           \
        DimensionedFieldType##SS##DimensionedFieldFunction;                    \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        DimensionedFieldType##SS##DimensionedFieldFunction,                    \
        0                                                                      \
    );                                                                         \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        DimensionedFieldType##DimensionedFieldFunction,                        \
        DimensionedFieldType##SS##DimensionedFieldFunction,                    \
        dictionary                                                             \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "DimensionedFieldFunction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
