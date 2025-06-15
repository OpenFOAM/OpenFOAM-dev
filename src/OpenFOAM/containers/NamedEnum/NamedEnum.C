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
Foam::NamedEnum<Enum, nEnum>::NamedEnum()
:
    HashTable<unsigned int>(2*nEnum)
{
    for (unsigned int ei = 0; ei < nEnum; ei++)
    {
        if (!names[ei] || names[ei][0] == '\0')
        {
            stringList goodNames(ei);

            for (unsigned int i = 0; i < ei; ++i)
            {
                goodNames[i] = names[i];
            }

            FatalErrorInFunction
                << "Illegal enumeration name at position " << ei << endl
                << "after entries " << goodNames << ".\n"
                << "Possibly your NamedEnum<Enum, nEnum>::names array"
                << " is not of size " << nEnum << endl
                << abort(FatalError);
        }
        insert(names[ei], ei);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Enum, unsigned int nEnum>
Enum Foam::NamedEnum<Enum, nEnum>::read(Istream& is) const
{
    const word name(is);

    HashTable<unsigned int>::const_iterator iter = find(name);

    if (iter == HashTable<unsigned int>::end())
    {
        FatalIOErrorInFunction(is)
            << name << " is not in enumeration: "
            << sortedToc() << exit(FatalIOError);
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
    Enum selection = Enum(0);
    unsigned int nSelections = 0;
    for (unsigned int ei = 0; ei < nEnum; ei++)
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
            << sortedToc() << exit(FatalIOError);
    }
    else if (nSelections > 1)
    {
        FatalIOErrorInFunction(dict)
            << "More than one option selected, please specify one of: "
            << sortedToc() << exit(FatalIOError);
    }

    return selection;
}


template<class Enum, unsigned int nEnum>
void Foam::NamedEnum<Enum, nEnum>::write(const Enum e, Ostream& os) const
{
    os  << operator[](e);
}


template<class Enum, unsigned int nEnum>
Foam::wordList Foam::NamedEnum<Enum, nEnum>::words()
{
    wordList lst(nEnum);

    label nElem = 0;
    for (unsigned int ei = 0; ei < nEnum; ei++)
    {
        if (names[ei] && names[ei][0])
        {
            lst[nElem++] = names[ei];
        }
    }

    lst.setSize(nElem);
    return lst;
}


template<class Enum, unsigned int nEnum>
const char* Foam::NamedEnum<Enum, nEnum>::operator[](const Enum e) const
{
    unsigned int ue = unsigned(e);

    if (ue < nEnum)
    {
        return names[ue];
    }
    else
    {
        FatalErrorInFunction
            << "names array index " << ue << " out of range 0-"
            << nEnum - 1
            << exit(FatalError);

        return names[0];
    }
}


// ************************************************************************* //
