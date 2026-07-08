/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2026 OpenFOAM Foundation
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
    Foam::Function1s::Antoine

Description
    Antoine equation, used to describe vapour pressure as a function of
    temperature

    \f[
        \log_{10} p = A - \frac{B}{C + T}
    \f]

    Additional variants are also provided which returns the natural logarithm
    of the pressure, and which evaluate the inverse of the function to obtain
    the temperature as a function of the pressure

SourceFiles
    Antoine.C

\*---------------------------------------------------------------------------*/

#ifndef Antoine_H
#define Antoine_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                       Class AntoineCoeffs Declaration
\*---------------------------------------------------------------------------*/

class AntoineCoeffs
{
protected:

    // Protected Data

        //- Natural logarithm of ten
        static const scalar ln10_;

        //- Coefficients
        const scalar A_, B_, C_;


public:

    // Constructors

        //- Construct from dictionary
        AntoineCoeffs
        (
            const word& typeName,
            const unitSets& units,
            const dictionary& dict
        );


    // Member Functions

        //- Write in dictionary format
        void write(Ostream& os, const unitSets& units) const;
};


/*---------------------------------------------------------------------------*\
                        Class AntoineFunction1 Declaration
\*---------------------------------------------------------------------------*/

template<class Function1Type>
class AntoineFunction1
:
    public FieldFunction1<scalar, Function1Type>,
    public AntoineCoeffs
{
public:

    // Constructors

        //- Construct from a name and dictionary
        inline AntoineFunction1
        (
            const word& name,
            const unitSets& units,
            const dictionary& dict
        );


    //- Destructor
    virtual inline ~AntoineFunction1();


    // Member Functions

        //- Write in dictionary format
        virtual inline void write(Ostream& os, const unitSets& units) const;
};


/*---------------------------------------------------------------------------*\
                           Class Antoine Declaration
\*---------------------------------------------------------------------------*/

class Antoine
:
    public AntoineFunction1<Antoine>
{
public:

    //- Runtime type information
    TypeName("Antoine");


    // Constructors

        //- Inherit Constructor
        using AntoineFunction1<Antoine>::AntoineFunction1;


    // Member Functions

        //- Return the value at a scalar coordinate
        virtual inline scalar value(const scalar T) const;

        //- Return the derivative at a scalar coordinate
        virtual inline scalar derivative(const scalar T) const;

        //- Return the integral between two scalar coordinates
        virtual inline scalar integral(const scalar T1, const scalar T2) const;
};


/*---------------------------------------------------------------------------*\
                           Class logAntoine Declaration
\*---------------------------------------------------------------------------*/

class logAntoine
:
    public AntoineFunction1<logAntoine>
{
public:

    //- Runtime type information
    TypeName("logAntoine");


    // Constructors

        //- Inherit Constructor
        using AntoineFunction1<logAntoine>::AntoineFunction1;


    // Member Functions

        //- Return the value at a scalar coordinate
        virtual inline scalar value(const scalar T) const;

        //- Return the derivative at a scalar coordinate
        virtual inline scalar derivative(const scalar T) const;

        //- Return the integral between two scalar coordinates
        virtual inline scalar integral(const scalar T1, const scalar T2) const;
};


/*---------------------------------------------------------------------------*\
                           Class inverseAntoine Declaration
\*---------------------------------------------------------------------------*/

class inverseAntoine
:
    public AntoineFunction1<inverseAntoine>
{
public:

    //- Runtime type information
    TypeName("inverseAntoine");


    // Constructors

        //- Inherit Constructor
        using AntoineFunction1<inverseAntoine>::AntoineFunction1;


    // Member Functions

        //- Return the value at a scalar coordinate
        virtual inline scalar value(const scalar p) const;

        //- Return the derivative at a scalar coordinate
        virtual inline scalar derivative(const scalar p) const;

        //- Return the integral between two scalar coordinates
        virtual inline scalar integral(const scalar p1, const scalar p2) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AntoineI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
