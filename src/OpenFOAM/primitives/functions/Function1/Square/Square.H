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
    Foam::Function1s::Square

Description
    Templated square-wave function with support for an offset level.

        \f[
            a square(f (x - x_0) ) + l
        \f]

    where

    \f$square(x)\f$ is the square-wave function in range \f$[-1, 1]\f$ with a
    mark/space ratio of \f$r\f$

    \vartable
        Symbol  | Description      | Data type       | Default
        a       | Amplitude        | Function1<Type> |
        f       | Frequency        | scalar          |
        x_0     | Phase            | scalar          | 0
        l       | Offset level     | Function1<Type> |
        r       | Mark/space ratio | scalar          | 1
    \endvartable

    Example for a scalar:
    \verbatim
        <name> square;
        <name>Coeffs
        {
            amplitude 2e-7;
            frequency 10;
            start     0;
            level     2e-6;
            markSpace 0.5;
        }
    \endverbatim

    Example for a vector:
    \verbatim
        <name> square;
        <name>Coeffs
        {
            amplitude (1 0.1 0);
            frequency 10;
            start     0;
            level     (10 1 0);
            markSpace 0.5;
        }
    \endverbatim

SourceFiles
    Square.C

\*---------------------------------------------------------------------------*/

#ifndef Square_H
#define Square_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                           Class Square Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Square
:
    public FieldFunction1<Type, Square<Type>>
{
    // Private Data

        //- Scalar amplitude of the square function
        autoPtr<Function1<Type>> amplitude_;

        //- Frequency of the square function
        scalar frequency_;

        //- Argument offset
        scalar start_;

        //- Level to which the square function is added
        autoPtr<Function1<Type>> level_;

        //- Mark/space ratio of the square function; the fraction of one period
        //  in which the function value is at the higher of the two possible
        //  values.
        scalar markSpace_;

        //- Is this function integrable?
        bool integrable_;


    // Private Member Functions

        //- Read the coefficients from the given dictionary
        void read(const dictionary& dict);


public:

    // Runtime type information
    TypeName("square");


    // Constructors

        //- Construct from name and dictionary
        Square
        (
            const word& name,
            const dictionary& dict
        );

        //- Copy constructor
        Square(const Square<Type>& se);


    //- Destructor
    virtual ~Square();


    // Member Functions

        //- Return value
        virtual inline Type value(const scalar x) const;

        //- Integrate between two values
        virtual inline Type integral(const scalar x1, const scalar x2) const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Square<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SquareI.H"

#ifdef NoRepository
    #include "Square.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
