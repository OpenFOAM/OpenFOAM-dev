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

#include "ListHashTable.H"
#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::ListHashTable<T, Key, Hash>::ListHashTable
(
    Istream& is,
    const label size
)
:
    ListHashTableCore(),
    keys_(ListHashTableCore::canonicalSize(size)),
    objects_(ListHashTableCore::canonicalSize(size)),
    nElmts_(0),
    endIter_(*this, keys_.size(), 0),
    endConstIter_(*this, keys_.size(), 0)
{
    if (size < 1)
    {
        FatalErrorInFunction
            << "Illegal size " << size << " for ListHashTable."
            << " Minimum size is 1" << abort(FatalError);
    }

    operator>>(is, *this);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::Ostream&
Foam::ListHashTable<T, Key, Hash>::printInfo(Ostream& os) const
{
    label used = 0;
    label maxChain = 0;
    unsigned avgChain = 0;

    // Find first non-empty entry
    forAll(keys_, hashIdx)
    {
        const label count = keys_[hashIdx].size();
        if (count)
        {
            ++used;
            avgChain += count;

            if (maxChain < count)
            {
                maxChain = count;
            }
        }
    }

    os  << "ListHashTable<T,Key,Hash>"
        << " elements:" << size() << " slots:" << used << "/" << keys_.size()
        << " chaining(avg/max):" << (used ? float(avgChain/used) : 0)
        << "/" << maxChain << endl;

    return os;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::Istream& Foam::operator>>(Istream& is, ListHashTable<T, Key, Hash>& L)
{
    is.fatalCheck("operator>>(Istream&, ListHashTable<T, Key, Hash>&)");

    // Anull list
    L.clear();

    is.fatalCheck("operator>>(Istream&, ListHashTable<T, Key, Hash>&)");

    token firstToken(is);

    is.fatalCheck
    (
        "operator>>(Istream&, ListHashTable<T, Key, Hash>&) : "
        "reading first token"
    );

    if (firstToken.isLabel())
    {
        label s = firstToken.labelToken();

        // Read beginning of contents
        char delimiter = is.readBeginList("ListHashTable<T, Key, Hash>");

        if (s)
        {
            if (2*s > L.keys_.size())
            {
                L.resize(2*s);
            }

            if (delimiter == token::BEGIN_LIST)
            {
                for (label i=0; i<s; i++)
                {
                    Key key;
                    is >> key;
                    L.insert(key, pTraits<T>(is));

                    is.fatalCheck
                    (
                        "operator>>(Istream&, ListHashTable<T, Key, Hash>&)"
                        " : reading entry"
                    );
                }
            }
            else
            {
                FatalIOErrorInFunction
                (
                    is
                )   << "incorrect first token, '(', found " << firstToken.info()
                    << exit(FatalIOError);
            }
        }

        // Read end of contents
        is.readEndList("ListHashTable");
    }
    else if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorInFunction
            (
                is
            )   << "incorrect first token, '(', found " << firstToken.info()
                << exit(FatalIOError);
        }

        token lastToken(is);
        while
        (
           !(
                lastToken.isPunctuation()
             && lastToken.pToken() == token::END_LIST
            )
        )
        {
            is.putBack(lastToken);

            Key key;
            is >> key;

            T element;
            is >> element;

            L.insert(key, element);

            is.fatalCheck
            (
                "operator>>(Istream&, ListHashTable<T, Key, Hash>&) : "
                "reading entry"
            );

            is >> lastToken;
        }
    }
    else
    {
        FatalIOErrorInFunction
        (
            is
        )   << "incorrect first token, expected <int> or '(', found "
            << firstToken.info()
            << exit(FatalIOError);
    }

    is.fatalCheck("operator>>(Istream&, ListHashTable<T, Key, Hash>&)");

    return is;
}


template<class T, class Key, class Hash>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ListHashTable<T, Key, Hash>& L)
{
    // Write size and start delimiter
    os << nl << L.size() << nl << token::BEGIN_LIST << nl;

    // Write contents
    for
    (
        typename ListHashTable<T, Key, Hash>::const_iterator iter = L.begin();
        iter != L.end();
        ++iter
    )
    {
        os << iter.key() << token::SPACE << iter() << nl;
    }

    // Write end delimiter
    os << token::END_LIST;

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, const ListHashTable&)");

    return os;
}


// ************************************************************************* //
