/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::complex

Description
    Extension to the c++ complex library type.

SourceFiles
    complexI.H
    complex.C

\*---------------------------------------------------------------------------*/

#ifndef complex_H
#define complex_H

#include "scalar.H"
#include "bool.H"
#include "word.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class complex;

inline scalar magSqr(const complex&);
inline complex sqr(const complex&);
inline scalar mag(const complex&);
inline const complex& max(const complex&, const complex&);
inline const complex& min(const complex&, const complex&);
inline complex limit(const complex&, const complex&);
inline const complex& sum(const complex&);
inline complex operator+(const complex&, const complex&);
inline complex operator-(const complex&);
inline complex operator-(const complex&, const complex&);
inline complex operator*(const complex&, const complex&);
inline complex operator/(const complex&, const complex&);
inline complex operator*(const scalar, const complex&);
inline complex operator*(const complex&, const scalar);
inline complex operator/(const complex&, const scalar);
inline complex operator/(const scalar, const complex&);
Istream& operator>>(Istream&, complex&);
Ostream& operator<<(Ostream&, const complex&);


/*---------------------------------------------------------------------------*\
                           Class complex Declaration
\*---------------------------------------------------------------------------*/

class complex
{
    // private data

        //- Real and imaginary parts of the complex number
        scalar re, im;

public:

    //- Component type
    typedef complex cmptType;


    // Static Data Members

        static const char* const typeName;

        static const complex zero;
        static const complex one;


    // Constructors

        //- Construct null
        inline complex();

        //- Construct given real and imaginary parts
        inline complex(const scalar Re, const scalar Im);

        //- Construct from Istream
        complex(Istream&);


    // Member Functions

           // Access

               inline scalar Re() const;
               inline scalar Im() const;

           // Edit

               inline scalar& Re();
               inline scalar& Im();

           // Operators

               inline complex conjugate() const;


    // Member Operators

        inline void operator+=(const complex&);
        inline void operator-=(const complex&);
        inline void operator*=(const complex&);
        inline void operator/=(const complex&);

        inline void operator=(const scalar);
        inline void operator+=(const scalar);
        inline void operator-=(const scalar);
        inline void operator*=(const scalar);
        inline void operator/=(const scalar);

        inline complex operator!() const;

        inline bool operator==(const complex&) const;
        inline bool operator!=(const complex&) const;


    // Friend Functions

        friend scalar magSqr(const complex& c);
        friend complex sqr(const complex& c);
        friend scalar mag(const complex& c);
        friend const complex& max(const complex&, const complex&);
        friend const complex& min(const complex&, const complex&);

        friend complex limit(const complex&, const complex&);

        friend const complex& sum(const complex&);


    // Friend operators

        friend complex operator+(const complex&, const complex&);
        friend complex operator-(const complex&);
        friend complex operator-(const complex&, const complex&);
        friend complex operator*(const complex&, const complex&);
        friend complex operator/(const complex&, const complex&);

        friend complex operator*(const scalar, const complex&);
        friend complex operator*(const complex&, const scalar);
        friend complex operator/(const complex&, const scalar);
        friend complex operator/(const scalar, const complex&);


    // IOstream Operators

        friend Istream& operator>>(Istream&, complex&);
        friend Ostream& operator<<(Ostream&, const complex&);

};


// * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * * //

//- Return a string representation of a complex
word name(const complex&);


//- Data associated with complex type are contiguous
template<>
inline bool contiguous<complex>() {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "complexI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
