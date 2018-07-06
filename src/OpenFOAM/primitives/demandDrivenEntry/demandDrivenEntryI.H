/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline void Foam::demandDrivenEntry<Type>::initialise() const
{
    if (!stored_)
    {
        dict_.lookup(keyword_) >> value_;
        stored_ = true;
    }
}


template<class Type>
inline const Type& Foam::demandDrivenEntry<Type>::value() const
{
    initialise();

    return value_;
}


template<class Type>
inline void Foam::demandDrivenEntry<Type>::setValue(const Type& value)
{
//    dict_.set<Type>(keyword_, value);
    value_ = value;
    stored_ = true;
}


template<class Type>
inline void Foam::demandDrivenEntry<Type>::reset()
{
    stored_ = false;
}


// ************************************************************************* //
