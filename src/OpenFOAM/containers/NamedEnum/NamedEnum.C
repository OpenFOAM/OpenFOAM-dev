/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "NamedEnum.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Enum, unsigned int nEnum>
Foam::NamedEnum<Enum, nEnum>::NamedEnum(const FixedList<word, nEnum>& names)
:
    FixedList<word, nEnum>(names),
    table_(2*nEnum)
{
    forAll(names, ei)
    {
        table_.insert(names[ei], ei);
    }
}


template<class Enum, unsigned int nEnum>
Foam::NamedEnum<Enum, nEnum>::NamedEnum(std::initializer_list<word> lst)
:
    NamedEnum(FixedList<word, nEnum>(lst))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Enum, unsigned int nEnum>
Enum Foam::NamedEnum<Enum, nEnum>::read(Istream& is) const
{
    const word name(is);

    HashTable<unsigned int>::const_iterator iter = table_.find(name);

    if (iter == HashTable<unsigned int>::end())
    {
        FatalIOErrorInFunction(is)
            << name << " is not in enumeration: "
            << *this << exit(FatalIOError);
    }

    return Enum(iter());
}


template<class Enum, unsigned int nEnum>
Enum Foam::NamedEnum<Enum, nEnum>::lookupOrDefault
(
    const word& name,
    const dictionary& dict,
    const Enum defaultValue
) const
{
    return dict.found(name)
      ? read(dict.lookup(name))
      : defaultValue;
}


template<class Enum, unsigned int nEnum>
Enum Foam::NamedEnum<Enum, nEnum>::select(const dictionary& dict) const
{
    const FixedList<word, nEnum>& names = *this;

    Enum selection = Enum(0);
    unsigned int nSelections = 0;
    forAll(names, ei)
    {
        if (dict.found(names[ei]))
        {
            selection = Enum(ei);
            nSelections++;
        }
    }

    if (nSelections == 0)
    {
        FatalIOErrorInFunction(dict)
            << "None of the options selected, please specify one of: "
            << names << exit(FatalIOError);
    }
    else if (nSelections > 1)
    {
        FatalIOErrorInFunction(dict)
            << "More than one option selected, please specify one of: "
            << names << exit(FatalIOError);
    }

    return selection;
}


template<class Enum, unsigned int nEnum>
void Foam::NamedEnum<Enum, nEnum>::write(const Enum e, Ostream& os) const
{
    os  << operator[](e);
}


// ************************************************************************* //
