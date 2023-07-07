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

InNamespace
    Foam

Description
    Functions for calculating the bounds and number of cells of a background
    mesh configured within a blockMeshDict file

\*---------------------------------------------------------------------------*/

#ifndef blockMeshFunctions_H
#define blockMeshFunctions_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline bool isEven(const label l)
{
    return (l % 2 == 0) ? true : false;
}


inline bool isEven(const Vector<label>& vl)
{
    return isEven(vl.x()) && isEven(vl.y()) && isEven(vl.z());
}


inline int order(const scalar s)
{
    return std::round(Foam::log10(mag(s)));
}


inline scalar roundingScale(const scalar s)
{
    return (s == 0) ? 0 : Foam::pow(10.0, order(s) - 1);
}


inline scalar roundDown(const scalar x, const scalar s)
{
    return floor(x/s)*s;
}


inline scalar roundUp(const scalar x, const scalar s)
{
    return ceil(x/s)*s;
}


inline vector roundDown(const vector& v, const scalar s)
{
    return vector
    (
        roundDown(v.x(), s),
        roundDown(v.y(), s),
        roundDown(v.z(), s)
    );
}


inline vector roundUp(const vector& v, const scalar s)
{
    return vector
    (
        roundUp(v.x(), s),
        roundUp(v.y(), s),
        roundUp(v.z(), s)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
