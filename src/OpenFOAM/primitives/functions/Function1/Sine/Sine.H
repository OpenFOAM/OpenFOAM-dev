/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    Foam::Function1Types::Sine

Description
    Templated sine function with support for an offset level.

        \f[
            a sin(2 \pi f (t - t_0)) s + l
        \f]

    where

    \vartable
        symbol  | Description       | Data type
        a       | Amplitude         | Function1<scalar>
        f       | Frequency [1/s]   | Function1<scalar>
        s       | Type scale factor | Function1<Type>
        l       | Type offset level | Function1<Type>
        t_0     | Start time [s]    | scalar
        t       | Time [s]          | scalar
    \endvartable

    Example for a scalar:
    \verbatim
        <entryName> sine;
        <entryName>Coeffs
        {
            frequency 10;
            amplitude 0.1;
            scale     2e-6;
            level     2e-6;
        }
    \endverbatim

    Example for a vector:
    \verbatim
        <entryName> sine;
        <entryName>Coeffs
        {
            frequency 10;
            amplitude 1;
            scale     (1 0.1 0);
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
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class Sine Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Sine
:
    public Function1<Type>
{
    // Private data

        //- Start-time for the sin function
        scalar t0_;

        //- Scalar amplitude of the sin function
        autoPtr<Function1<scalar>> amplitude_;

        //- Frequency of the sin function
        autoPtr<Function1<scalar>> frequency_;

        //- Scaling factor of the sin function
        autoPtr<Function1<Type>> scale_;

        //- Level to which the sin function is added
        autoPtr<Function1<Type>> level_;


    // Private Member Functions

        //- Read the coefficients from the given dictionary
        void read(const dictionary& coeffs);

        //- Disallow default bitwise assignment
        void operator=(const Sine<Type>&);


public:

    // Runtime type information
    TypeName("sine");


    // Constructors

        //- Construct from entry name and dictionary
        Sine
        (
            const word& entryName,
            const dictionary& dict
        );

        //- Copy constructor
        Sine(const Sine<Type>& se);


    //- Destructor
    virtual ~Sine();


    // Member Functions

        //- Return value for time t
        virtual inline Type value(const scalar t) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SineI.H"

#ifdef NoRepository
    #include "Sine.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
