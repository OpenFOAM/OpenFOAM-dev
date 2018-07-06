/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Enum, unsigned int nEnum>
Foam::NamedEnum<Enum, nEnum>::NamedEnum()
:
    HashTable<unsigned int>(2*nEnum)
{
    for (unsigned int enumI = 0; enumI < nEnum; ++enumI)
    {
        if (!names[enumI] || names[enumI][0] == '\0')
        {
            stringList goodNames(enumI);

            for (unsigned int i = 0; i < enumI; ++i)
            {
                goodNames[i] = names[i];
            }

            FatalErrorInFunction
                << "Illegal enumeration name at position " << enumI << endl
                << "after entries " << goodNames << ".\n"
                << "Possibly your NamedEnum<Enum, nEnum>::names array"
                << " is not of size " << nEnum << endl
                << abort(FatalError);
        }
        insert(names[enumI], enumI);
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
void Foam::NamedEnum<Enum, nEnum>::write(const Enum e, Ostream& os) const
{
    os  << operator[](e);
}


template<class Enum, unsigned int nEnum>
Foam::stringList Foam::NamedEnum<Enum, nEnum>::strings()
{
    stringList lst(nEnum);

    label nElem = 0;
    for (unsigned int enumI = 0; enumI < nEnum; ++enumI)
    {
        if (names[enumI] && names[enumI][0])
        {
            lst[nElem++] = names[enumI];
        }
    }

    lst.setSize(nElem);
    return lst;
}


template<class Enum, unsigned int nEnum>
Foam::wordList Foam::NamedEnum<Enum, nEnum>::words()
{
    wordList lst(nEnum);

    label nElem = 0;
    for (unsigned int enumI = 0; enumI < nEnum; ++enumI)
    {
        if (names[enumI] && names[enumI][0])
        {
            lst[nElem++] = names[enumI];
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
