/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::Function1s::Dimensioned

Description
    Function1 with dimensions.

    The dimensions are specified in the constructing code, in much the same way
    as for a dimensioned constant, except that both argument (x) and value
    dimensions must be specified.

    This can be used identically to a standard Function1. In addition, it also
    permits specification of the dimensions. This allows the dimensions to be
    checked. It also permits unit conversions.

Usage
    <name>
    {
        type        table;

        // Dimensions
        xDimensions [CAD];  // Optional. Argument dimensions/units.
                            // Here, this specifies coordinates are in
                            // crank angle degrees (available if using
                            // engine time).
        dimensions  [mm];   // Optional. Value dimensions/units.
                            // Here, values are in mm.

        // Tabulated data
        values
        (
            (0 0)
            (60 12)         // <-- i.e., 12 mm at 60 degrees
            (180 20)
            (240 8)
            (360 0)
        );
        outOfBounds repeat;
    }

SourceFiles
    DimensionedFunction1.C
    DimensionedFunction1I.H

\*---------------------------------------------------------------------------*/

#ifndef DimensionedFunction1_H
#define DimensionedFunction1_H

#include "Function1.H"
#include "dimensionedType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                            Class Dimensioned Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Dimensioned
:
    public FieldFunction1<Type, Dimensioned<Type>>
{
    // Private Data

        //- Argument dimensions
        dimensionedScalar xUnits_;

        //- Value dimensions
        dimensionedScalar units_;

        //- Value function
        autoPtr<Function1<Type>> value_;


    // Private Member Functions

        //- Read the coefficients from the given dictionary. Note that this is
        //  the outer dictionary, not the coefficients dictionary, as is the
        //  case for most/all other Function1-s.
        void read(const word& name, const dictionary& dict);


public:

    // Runtime type information
    virtual const word& type() const
    {
        return value_->type();
    }


    // Constructors

        //- Construct from name and dictionary. Note that this is the outer
        //  dictionary, not the coefficients dictionary, as is the case for
        //  most/all other Function1-s.
        Dimensioned
        (
            const word& name,
            const dimensionSet& xDimensions,
            const dimensionSet& dimensions,
            const dictionary& dict
        );

        //- Copy constructor
        Dimensioned(const Dimensioned<Type>& se);


    //- Destructor
    virtual ~Dimensioned();


    // Member Functions

        //- Return a value
        virtual inline Type value(const scalar x) const;

        //- Integrate between two values
        virtual inline Type integral(const scalar x1, const scalar x2) const;

        //- Inherit field value evaluation
        using FieldFunction1<Type, Dimensioned<Type>>::value;

        //- Inherit field integral evaluation
        using FieldFunction1<Type, Dimensioned<Type>>::integral;

        //- Return a dimensioned value
        virtual inline dimensioned<Type> value
        (
            const dimensionedScalar& x
        ) const;

        //- Integrate between two dimensioned values
        virtual inline dimensioned<Type> integral
        (
            const dimensionedScalar& x1,
            const dimensionedScalar& x2
        ) const;

        //- Is this function guaranteed to be constant?
        virtual inline bool constant() const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Dimensioned<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DimensionedFunction1I.H"

#ifdef NoRepository
    #include "DimensionedFunction1.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
