/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "UniformTable2.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2s::UniformTable<Type>::UniformTable
(
    const word& name,
    const dictionary& dict
)
:
    FieldFunction2<Type, UniformTable<Type>>(name),
    low_(dict.lookup("low")),
    high_(dict.lookup("high")),
    values_(dict.lookup("values"))
{
    if (values_.m() < 2 || values_.n() < 2)
    {
        FatalErrorInFunction
            << "Table " << nl
            << "    " << this->name_ << nl
            << "    has less than 2 entries in one or both dimensions."
            << exit(FatalError);
    }
    else
    {
        deltax_ = (high_.first() - low_.first())/(values_.m() - 1);
        deltay_ = (high_.second() - low_.second())/(values_.n() - 1);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline void Foam::Function2s::UniformTable<Type>::checkRange
(
    scalar x,
    scalar ndx,
    label ix,
    scalar y,
    scalar ndy,
    label iy
) const
{
    if (ndx < 0 || ix > values_.m() - 2)
    {
        FatalErrorInFunction
            << "x " << x << " out of range "
            << low_.first() << " to " << high_.first() << nl
            << "    of table " << this->name_
            << exit(FatalError);
    }

    if (ndy < 0 || iy > values_.n() - 2)
    {
        FatalErrorInFunction
            << "y " << y << " out of range "
            << low_.second() << " to " << high_.second() << nl
            << "    of table " << this->name_
            << exit(FatalError);
    }
}


template<class Type>
Type Foam::Function2s::UniformTable<Type>::value
(
    scalar x,
    scalar y
) const
{
    const scalar ndx = (x - low_.first())/deltax_;
    const label ix = ndx;

    const scalar ndy = (y - low_.second())/deltay_;
    const label iy = ndy;

    checkRange(x, ndx, ix, y, ndy, iy);

    const scalar xi = low_.first() + ix*deltax_;
    const scalar lambdax = (x - xi)/deltax_;

    // Interpolate the values at yi wrt x
    const Type fxi =
        values_(ix, iy)
      + lambdax*(values_(ix + 1, iy) - values_(ix, iy));

    // Interpolate the values at yi+1 wrt x
    const Type fxix1 =
        values_(ix, iy + 1)
      + lambdax*(values_(ix + 1, iy + 1) - values_(ix, iy + 1));

    const scalar yi = low_.second() + iy*deltay_;
    const scalar lambday = (y - yi)/deltay_;

    // Interpolate wrt y
    return fxi + lambday*(fxix1 - fxi);
}


template<class Type>
Type Foam::Function2s::UniformTable<Type>::
dfdp
(
    scalar p,
    scalar T
) const
{
    const scalar ndp = (p - low_.first())/deltax_;
    const label ip = ndp;

    const scalar ndT = (T - low_.second())/deltay_;
    const label iT = ndT;

    checkRange(p, ndp, ip, T, ndT, iT);

    const Type dfdpi =
        (values_(ip + 1, iT) - values_(ip, iT))/deltax_;
    const Type dfdpip1 =
        (values_(ip + 1, iT + 1) - values_(ip, iT + 1))/deltax_;

    const scalar Ti = low_.second() + iT*deltay_;
    const scalar lambdaT = (T - Ti)/deltay_;

    // Interpolate wrt T
    return dfdpi + lambdaT*(dfdpip1 - dfdpi);
}


template<class Type>
Type Foam::Function2s::UniformTable<Type>::
dfdT
(
    scalar p,
    scalar T
) const
{
    const scalar ndp = (p - low_.first())/deltax_;
    const label ip = ndp;

    const scalar ndT = (T - low_.second())/deltay_;
    const label iT = ndT;

    checkRange(p, ndp, ip, T, ndT, iT);

    const Type dfdTi =
        (values_(ip, iT + 1) - values_(ip, iT))/deltay_;
    const Type dfdTip1 =
        (values_(ip + 1, iT + 1) - values_(ip + 1, iT))/deltay_;

    const scalar pi = low_.first() + ip*deltax_;
    const scalar lambdap = (p - pi)/deltax_;

    // Interpolate wrt p
    return dfdTi + lambdap*(dfdTip1 - dfdTi);
}


template<class Type>
void Foam::Function2s::UniformTable<Type>::write(Ostream& os) const
{
    writeEntry(os, "low", low_);
    writeEntry(os, "high", high_);
    writeEntry(os, "values", values_);
}


// ************************************************************************* //
