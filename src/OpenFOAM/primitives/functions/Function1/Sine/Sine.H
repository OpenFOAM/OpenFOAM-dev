/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
        a       | Amplitude
        f       | Frequency [1/s]
        s       | Type scale factor
        l       | Type offset level
        t_0     | Start time [s]
        t       | Time [s]
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

// Forward declaration of friend functions and operators
namespace Function1Types
{
    template<class Type> class Sine;
};

template<class Type>
Ostream& operator<<(Ostream&, const Function1Types::Sine<Type>&);

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
        scalar amplitude_;

        //- Frequency of the sin function
        scalar frequency_;

        //- Scaling factor of the sin function
        Type scale_;

        //- Level to which the sin function is added
        Type level_;


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
            const dictionary& dict,
            const word& ext = "Coeffs"
        );

        //- Copy constructor
        Sine(const Sine<Type>& se);

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return tmp<Function1<Type>>(new Sine<Type>(*this));
        }


    //- Destructor
    virtual ~Sine();


    // Member Functions

        //- Return value for time t
        Type value(const scalar t) const;

        //- Integrate between the two time values t1 and t2
        Type integrate(const scalar t1, const scalar t2) const;


    // I/O

        //- Ostream Operator
        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const Sine<Type>& cnst
        );

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Sine.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
