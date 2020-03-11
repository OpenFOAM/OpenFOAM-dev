/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "UList.H"
#include "Ostream.H"
#include "token.H"
#include "SLList.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * IOstream Functions  * * * * * * * * * * * * //

template<class ListType>
void Foam::writeListEntry(Ostream& os, const ListType& l)
{
    if
    (
        token::compound::isCompound
        (
            "List<"
          + word(pTraits<typename ListType::value_type>::typeName) + '>'
        )
    )
    {
        os << word
        (
            "List<"
          + word(pTraits<typename ListType::value_type>::typeName) + '>'
        ) << " ";
    }

    os << l;
}


template<class ListType>
void Foam::writeListEntries(Ostream& os, const ListType& l)
{
    // Write size and start delimiter
    os << nl << l.size() << nl << token::BEGIN_LIST;

    // Write contents
    forAll(l, i)
    {
        writeEntry(os, l[i]);
        os << nl;
    }

    // Write end delimiter
    os << nl << token::END_LIST << nl;
}


template<class ListType>
void Foam::writeListEntries(Ostream& os, const word& keyword, const ListType& l)
{
    writeKeyword(os, keyword);
    writeListEntries(os, l);
    os << token::END_STATEMENT << endl;
}


template<class T>
void Foam::writeEntry(Ostream& os, const UList<T>& l)
{
    writeListEntry(os, l);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class T>
Foam::Ostream& Foam::operator<<(Foam::Ostream& os, const Foam::UList<T>& L)
{
    // Write list contents depending on data format
    if (os.format() == IOstream::ASCII || !contiguous<T>())
    {
        bool uniform = false;

        if (L.size() > 1 && contiguous<T>())
        {
            uniform = true;

            forAll(L, i)
            {
                if (L[i] != L[0])
                {
                    uniform = false;
                    break;
                }
            }
        }

        if (uniform)
        {
            // Write size and start delimiter
            os << L.size() << token::BEGIN_BLOCK;

            // Write contents
            os << L[0];

            // Write end delimiter
            os << token::END_BLOCK;
        }
        else if (L.size() <= 1 || (L.size() < 11 && contiguous<T>()))
        {
            // Write size and start delimiter
            os << L.size() << token::BEGIN_LIST;

            // Write contents
            forAll(L, i)
            {
                if (i > 0) os << token::SPACE;
                os << L[i];
            }

            // Write end delimiter
            os << token::END_LIST;
        }
        else
        {
            // Write size and start delimiter
            os << nl << L.size() << nl << token::BEGIN_LIST;

            // Write contents
            forAll(L, i)
            {
                os << nl << L[i];
            }

            // Write end delimiter
            os << nl << token::END_LIST << nl;
        }
    }
    else
    {
        os << nl << L.size() << nl;
        if (L.size())
        {
            os.write(reinterpret_cast<const char*>(L.v_), L.byteSize());
        }
    }

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, const UList&)");

    return os;
}


template<class T>
Foam::Istream& Foam::operator>>(Istream& is, UList<T>& L)
{
    is.fatalCheck("operator>>(Istream&, UList<T>&)");

    token firstToken(is);

    is.fatalCheck("operator>>(Istream&, UList<T>&) : reading first token");

    if (firstToken.isCompound())
    {
        List<T> elems;
        elems.transfer
        (
            dynamicCast<token::Compound<List<T>>>
            (
                firstToken.transferCompoundToken(is)
            )
        );
        // Check list length
        label s = elems.size();

        if (s != L.size())
        {
            FatalIOErrorInFunction(is)
                << "incorrect length for UList. Read " << s
                << " expected " << L.size()
                << exit(FatalIOError);
        }
        for (label i=0; i<s; i++)
        {
            L[i] = elems[i];
        }
    }
    else if (firstToken.isLabel())
    {
        label s = firstToken.labelToken();

        // Set list length to that read
        if (s != L.size())
        {
            FatalIOErrorInFunction(is)
                << "incorrect length for UList. Read " << s
                << " expected " << L.size()
                << exit(FatalIOError);
        }

        // Read list contents depending on data format

        if (is.format() == IOstream::ASCII || !contiguous<T>())
        {
            // Read beginning of contents
            char delimiter = is.readBeginList("List");

            if (s)
            {
                if (delimiter == token::BEGIN_LIST)
                {
                    for (label i=0; i<s; i++)
                    {
                        is >> L[i];

                        is.fatalCheck
                        (
                            "operator>>(Istream&, UList<T>&) : reading entry"
                        );
                    }
                }
                else
                {
                    T element;
                    is >> element;

                    is.fatalCheck
                    (
                        "operator>>(Istream&, UList<T>&) : "
                        "reading the single entry"
                    );

                    for (label i=0; i<s; i++)
                    {
                        L[i] = element;
                    }
                }
            }

            // Read end of contents
            is.readEndList("List");
        }
        else
        {
            if (s)
            {
                is.read(reinterpret_cast<char*>(L.data()), s*sizeof(T));

                is.fatalCheck
                (
                    "operator>>(Istream&, UList<T>&) : reading the binary block"
                );
            }
        }
    }
    else if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorInFunction(is)
                << "incorrect first token, expected '(', found "
                << firstToken.info()
                << exit(FatalIOError);
        }

        // Putback the opening bracket
        is.putBack(firstToken);

        // Now read as a singly-linked list
        SLList<T> sll(is);

        if (sll.size() != L.size())
        {
            FatalIOErrorInFunction(is)
                << "incorrect length for UList. Read " << sll.size()
                << " expected " << L.size()
                << exit(FatalIOError);
        }

        // Convert the singly-linked list to this list
        label i = 0;
        for
        (
            typename SLList<T>::const_iterator iter = sll.begin();
            iter != sll.end();
            ++iter
        )
        {
            L[i] = iter();
        }

    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int> or '(', found "
            << firstToken.info()
            << exit(FatalIOError);
    }

    return is;
}


// ************************************************************************* //
