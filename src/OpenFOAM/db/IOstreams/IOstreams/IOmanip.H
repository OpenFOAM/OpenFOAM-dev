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

InNamespace
    Foam::IOmanip

Description
    Istream and Ostream manipulators taking arguments.

\*---------------------------------------------------------------------------*/

#ifndef IOmanip_H
#define IOmanip_H

#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T> class Smanip;
template<class T> class Imanip;
template<class T> class Omanip;

template<class T>
inline Istream& operator>>(Istream& is, const Smanip<T>& m);

template<class T>
inline Ostream& operator<<(Ostream& os, const Smanip<T>& m);

template<class T>
inline Istream& operator>>(Istream& is, const Imanip<T>& m);

template<class T>
inline Ostream& operator<<(Ostream& os, const Omanip<T>& m);


/*---------------------------------------------------------------------------*\
                        Class Smanip Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class Smanip
{
    T (IOstream::*_fPtr)(const T);
    T _i;

public:

    Smanip(T (IOstream::*fPtr)(const T), const T i)
    :
        _fPtr(fPtr),
        _i(i)
    {}

    friend Istream& operator>> <T>(Istream& is, const Smanip<T>& m);
    friend Ostream& operator<< <T>(Ostream& os, const Smanip<T>& m);
};


template<class T>
inline Istream& operator>>(Istream& is, const Smanip<T>& m)
{
    (is.*m._fPtr)(m._i);
    return is;
}


template<class T>
inline Ostream& operator<<(Ostream& os, const Smanip<T>& m)
{
    (os.*m._fPtr)(m._i);
    return os;
}


/*---------------------------------------------------------------------------*\
                        Class Imanip Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class Imanip
{
    T (Istream::*_fPtr)(const T);
    T _i;

public:

    Imanip(T (Istream::*fPtr)(const T), const T i)
    :
        _fPtr(fPtr),
        _i(i)
    {}

    friend Istream& operator>> <T>(Istream& is, const Imanip<T>& m);
};


template<class T>
inline Istream& operator>>(Istream& is, const Imanip<T>& m)
{
    (is.*m._fPtr)(m._i);
    return is;
}


/*---------------------------------------------------------------------------*\
                        Class Omanip Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class Omanip
{
    T (Ostream::*_fPtr)(const T);
    T _i;

public:

    Omanip(T (Ostream::*fPtr)(const T), const T i)
    :
        _fPtr(fPtr),
        _i(i)
    {}

    friend Ostream& operator<< <T>(Ostream& os, const Omanip<T>& m);
};


template<class T>
inline Ostream& operator<<(Ostream& os, const Omanip<T>& m)
{
    (os.*m._fPtr)(m._i);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Smanip<ios_base::fmtflags> setf
(
    const ios_base::fmtflags flags
)
{
    return Smanip<ios_base::fmtflags>(&IOstream::setf, flags);
}


inline Omanip<IOstream::streamFormat> setformat
(
    const IOstream::streamFormat fmt
)
{
    return Omanip<IOstream::streamFormat>(&IOstream::format, fmt);
}


inline Omanip<IOstream::versionNumber> setversion
(
    const IOstream::versionNumber ver
)
{
    return Omanip<IOstream::versionNumber>(&IOstream::version, ver);
}


inline Omanip<IOstream::compressionType> setcompression
(
    const IOstream::compressionType cmp
)
{
    return Omanip<IOstream::compressionType>(&IOstream::compression, cmp);
}


inline Omanip<int> setw(const int i)
{
    return Omanip<int>(&Ostream::width, i);
}


inline Omanip<int> setprecision(const int i)
{
    return Omanip<int>(&Ostream::precision, i);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
