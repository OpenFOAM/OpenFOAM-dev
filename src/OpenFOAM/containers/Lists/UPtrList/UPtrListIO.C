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

\*---------------------------------------------------------------------------*/

#include "UPtrList.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::UPtrList<T>::writeEntry(Ostream& os) const
{
    if
    (
        size()
     && token::compound::isCompound
        (
            "List<" + word(pTraits<T>::typeName) + '>'
        )
    )
    {
        os  << word("List<" + word(pTraits<T>::typeName) + '>') << " ";
    }

    os << *this;
}


template<class T>
void Foam::UPtrList<T>::writeEntry(const word& keyword, Ostream& os) const
{
    os.writeKeyword(keyword);
    writeEntry(os);
    os << token::END_STATEMENT << endl;
}


template<class T>
void Foam::UPtrList<T>::writeEntryList(Ostream& os) const
{
    // Write size and start delimiter
    os << nl << size() << nl << token::BEGIN_LIST;

    // Write contents
    forAll(*this, i)
    {
        this->operator[](i).writeEntry(os);
        os << nl;
    }

    // Write end delimiter
    os << nl << token::END_LIST << nl;
}


template<class T>
void Foam::UPtrList<T>::writeEntryList(const word& keyword, Ostream& os) const
{
    os.writeKeyword(keyword);
    writeEntryList(os);
    os << token::END_STATEMENT << endl;
}


// * * * * * * * * * * * * * * * Ostream Operators * * * * * * * * * * * * * //

template<class T>
Foam::Ostream& Foam::operator<<(Ostream& os, const UPtrList<T>& L)
{
    // Write size and start delimiter
    os  << nl << indent << L.size() << nl
        << indent << token::BEGIN_LIST << incrIndent;

    // Write contents
    forAll(L, i)
    {
        os << nl << L[i];
    }

    // Write end delimiter
    os << nl << decrIndent << indent << token::END_LIST << nl;

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, const UPtrList&)");

    return os;
}


// ************************************************************************* //
