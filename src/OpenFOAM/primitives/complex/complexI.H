/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline complex::complex()
{}


inline complex::complex(const scalar Re, const scalar Im)
:
    re(Re),
    im(Im)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline scalar complex::Re() const
{
    return re;
}


inline scalar complex::Im() const
{
    return im;
}


inline scalar& complex::Re()
{
    return re;
}


inline scalar& complex::Im()
{
    return im;
}


inline complex complex::conjugate() const
{
    return complex(re, -im);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline const complex& complex::operator=(const complex& c)
{
    re = c.re;
    im = c.im;
    return *this;
}


inline void complex::operator+=(const complex& c)
{
    re += c.re;
    im += c.im;
}


inline void complex::operator-=(const complex& c)
{
    re -= c.re;
    im -= c.im;
}


inline void complex::operator*=(const complex& c)
{
    *this = (*this)*c;
}


inline void complex::operator/=(const complex& c)
{
    *this = *this/c;
}


inline const complex& complex::operator=(const scalar s)
{
    re = s;
    im = 0.0;
    return *this;
}


inline void complex::operator+=(const scalar s)
{
    re += s;
}


inline void complex::operator-=(const scalar s)
{
    re -= s;
}


inline void complex::operator*=(const scalar s)
{
    re *= s;
    im *= s;
}


inline void complex::operator/=(const scalar s)
{
    re /= s;
    im /= s;
}


inline complex complex::operator!() const
{
    return conjugate();
}


inline bool complex::operator==(const complex& c) const
{
    return (equal(re, c.re) && equal(im, c.im));
}


inline bool complex::operator!=(const complex& c) const
{
    return !operator==(c);
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


inline scalar magSqr(const complex& c)
{
    return (c.re*c.re + c.im*c.im);
}


inline complex sqr(const complex& c)
{
    return c * c;
}


inline scalar mag(const complex& c)
{
    return sqrt(magSqr(c));
}


inline const complex& max(const complex& c1, const complex& c2)
{
    if (mag(c1) > mag(c2))
    {
        return c1;
    }
    else
    {
        return c2;
    }
}


inline const complex& min(const complex& c1, const complex& c2)
{
    if (mag(c1) < mag(c2))
    {
        return c1;
    }
    else
    {
        return c2;
    }
}


inline complex limit(const complex& c1, const complex& c2)
{
    return complex(limit(c1.re, c2.re), limit(c1.im, c2.im));
}


inline const complex& sum(const complex& c)
{
    return c;
}


template<class Cmpt>
class Tensor;

inline complex transform(const Tensor<scalar>&, const complex c)
{
    return c;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline complex operator+(const complex& c1, const complex& c2)
{
    return complex
    (
        c1.re + c2.re,
        c1.im + c2.im
    );
}


inline complex operator-(const complex& c)
{
    return complex
    (
        -c.re,
        -c.im
    );
}


inline complex operator-(const complex& c1, const complex& c2)
{
    return complex
    (
        c1.re - c2.re,
        c1.im - c2.im
    );
}


inline complex operator*(const complex& c1, const complex& c2)
{
    return complex
    (
        c1.re*c2.re - c1.im*c2.im,
        c1.im*c2.re + c1.re*c2.im
    );
}


inline complex operator/(const complex& c1, const complex& c2)
{
    scalar sqrC2 = magSqr(c2);

    return complex
    (
        (c1.re*c2.re + c1.im*c2.im)/sqrC2,
        (c1.im*c2.re - c1.re*c2.im)/sqrC2
    );
}


inline complex operator*(const scalar s, const complex& c)
{
    return complex(s*c.re, s*c.im);
}


inline complex operator*(const complex& c, const scalar s)
{
    return complex(s*c.re, s*c.im);
}


inline complex operator/(const complex& c, const scalar s)
{
    return complex(c.re/s, c.im/s);
}


inline complex operator/(const scalar s, const complex& c)
{
    return complex(s/c.re, s/c.im);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
