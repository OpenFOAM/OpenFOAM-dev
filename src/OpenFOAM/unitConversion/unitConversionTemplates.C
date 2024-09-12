/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "unitConversion.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
T Foam::unitConversion::toStandard(const T& t) const
{
    return standard() ? t : t*multiplier_;
}


template<class T>
Foam::Pair<T> Foam::unitConversion::toStandard(const Pair<T>& p) const
{
    return Pair<T>(toStandard(p.first(), p.second()));
}


template<class T>
Foam::List<T> Foam::unitConversion::toStandard(const List<T>& l) const
{
    List<T> result(l.size());
    forAll(l, i)
    {
        result[i] = toStandard(l[i]);
    }
    return result;
}


template<class T>
Foam::tmp<Foam::Field<T>>
Foam::unitConversion::toStandard(const Field<T>& f) const
{
    return standard() ? tmp<Field<T>>(f) : f*multiplier_;
}


template<class T>
Foam::tmp<Foam::Field<T>>
Foam::unitConversion::toStandard(const tmp<Field<T>>& tf) const
{
    return standard() ? tf : tf*multiplier_;
}


template<class T>
void Foam::unitConversion::makeStandard(T& t) const
{
    if (!standard())
    {
        t *= multiplier_;
    }
}


template<class T>
void Foam::unitConversion::makeStandard(Pair<T>& p) const
{
    if (!standard())
    {
        makeStandard(p.first());
        makeStandard(p.second());
    }
}


template<class T>
void Foam::unitConversion::makeStandard(List<T>& l) const
{
    if (!standard())
    {
        forAll(l, i)
        {
            makeStandard(l[i]);
        }
    }
}


template<class T>
T Foam::unitConversion::toUser(const T& t) const
{
    return standard() ? t : t/multiplier_;
}


template<class T>
Foam::Pair<T> Foam::unitConversion::toUser(const Pair<T>& p) const
{
    return Pair<T>(toUser(p.first()), toUser(p.second()));
}


template<class T>
Foam::List<T> Foam::unitConversion::toUser(const List<T>& l) const
{
    List<T> result(l.size());
    forAll(l, i)
    {
        result[i] = toUser(l[i]);
    }
    return result;
}


template<class T>
Foam::tmp<Foam::Field<T>>
Foam::unitConversion::toUser(const Field<T>& f) const
{
    return standard() ? tmp<Field<T>>(f) : f/multiplier_;
}


template<class T>
Foam::tmp<Foam::Field<T>>
Foam::unitConversion::toUser(const tmp<Field<T>>& tf) const
{
    return standard() ? tf : tf/multiplier_;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::writeEntry
(
    Ostream& os,
    const unitConversion& defaultUnits,
    const Type& t
)
{
    return writeEntry(os, defaultUnits.toUser(t));
}


// ************************************************************************* //
