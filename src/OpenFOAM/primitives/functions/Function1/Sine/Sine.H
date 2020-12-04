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
    Foam::Function1s::Sine

Description
    Templated sine function with support for an offset level.

        \f[
            a sin(2 \pi f (x - x_0)) + l
        \f]

    where

    \vartable
        Symbol  | Description  | Data type       | Default
        a       | Amplitude    | Function1<Type> |
        f       | Frequency    | scalar          |
        x_0     | Start        | scalar          | 0
        l       | Offset level | Function1<Type> |
    \endvartable

    Example for a scalar:
    \verbatim
        <name> sine;
        <name>Coeffs
        {
            amplitude 2e-7;
            frequency 10;
            start     0;
            level     2e-6;
        }
    \endverbatim

    Example for a vector:
    \verbatim
        <name> sine;
        <name>Coeffs
        {
            amplitude (1 0.1 0);
            frequency 10;
            start     0;
            level     (10 1 0);
        }
    \endverbatim

SourceFiles
    Sine.C

\*---------------------------------------------------------------------------*/

#ifndef Sine_H
#define Sine_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                            Class Sine Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Sine
:
    public FieldFunction1<Type, Sine<Type>>
{
    // Private Data

        //- Amplitude of the sine function
        autoPtr<Function1<Type>> amplitude_;

        //- Frequency of the sine function
        scalar frequency_;

        //- Argument offset
        scalar start_;

        //- Level to which the sine function is added
        autoPtr<Function1<Type>> level_;

        //- Is this function integrable?
        bool integrable_;


    // Private Member Functions

        //- Read the coefficients from the given dictionary
        void read(const dictionary& dict);


public:

    // Runtime type information
    TypeName("sine");


    // Constructors

        //- Construct from name and dictionary
        Sine
        (
            const word& name,
            const dictionary& dict
        );

        //- Copy constructor
        Sine(const Sine<Type>& se);


    //- Destructor
    virtual ~Sine();


    // Member Functions

        //- Return value
        virtual inline Type value(const scalar x) const;

        //- Integrate between two values
        virtual inline Type integral(const scalar x1, const scalar x2) const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Sine<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SineI.H"

#ifdef NoRepository
    #include "Sine.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
