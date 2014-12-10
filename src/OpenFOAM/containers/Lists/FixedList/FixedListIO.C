/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "FixedList.H"
#include "Istream.H"
#include "Ostream.H"
#include "token.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, unsigned Size>
Foam::FixedList<T, Size>::FixedList(Istream& is)
{
    operator>>(is, *this);
}


template<class T, unsigned Size>
Foam::Istream& Foam::operator>>(Foam::Istream& is, FixedList<T, Size>& L)
{
    is.fatalCheck("operator>>(Istream&, FixedList<T, Size>&)");

    if (is.format() == IOstream::ASCII || !contiguous<T>())
    {
        token firstToken(is);

        is.fatalCheck
        (
            "operator>>(Istream&, FixedList<T, Size>&) : reading first token"
        );

        if (firstToken.isCompound())
        {
            L = dynamicCast<token::Compound<List<T> > >
            (
                firstToken.transferCompoundToken(is)
            );
        }
        else if (firstToken.isLabel())
        {
            label s = firstToken.labelToken();

            // Set list length to that read
            L.checkSize(s);
        }
        else if (!firstToken.isPunctuation())
        {
            FatalIOErrorIn("operator>>(Istream&, FixedList<T, Size>&)", is)
                << "incorrect first token, expected <label> "
                   "or '(' or '{', found "
                << firstToken.info()
                << exit(FatalIOError);
        }
        else
        {
            // Putback the opening bracket
            is.putBack(firstToken);
        }

        // Read beginning of contents
        char delimiter = is.readBeginList("FixedList");

        if (delimiter == token::BEGIN_LIST)
        {
            for (register unsigned i=0; i<Size; i++)
            {
                is >> L[i];

                is.fatalCheck
                (
                    "operator>>(Istream&, FixedList<T, Size>&) : "
                    "reading entry"
                );
            }
        }
        else
        {
            T element;
            is >> element;

            is.fatalCheck
            (
                "operator>>(Istream&, FixedList<T, Size>&) : "
                "reading the single entry"
            );

            for (register unsigned i=0; i<Size; i++)
            {
                L[i] = element;
            }
        }

        // Read end of contents
        is.readEndList("FixedList");
    }
    else
    {
        is.read(reinterpret_cast<char*>(L.data()), Size*sizeof(T));

        is.fatalCheck
        (
            "operator>>(Istream&, FixedList<T, Size>&) : "
            "reading the binary block"
        );
    }

    return is;
}


// * * * * * * * * * * * * * * * Ostream Operator *  * * * * * * * * * * * * //

template<class T, unsigned Size>
void Foam::FixedList<T, Size>::writeEntry(Ostream& os) const
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


template<class T, unsigned Size>
void Foam::FixedList<T, Size>::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    os.writeKeyword(keyword);
    writeEntry(os);
    os << token::END_STATEMENT << endl;
}


template<class T, unsigned Size>
Foam::Ostream& Foam::operator<<(Ostream& os, const FixedList<T, Size>& L)
{
    // Write list contents depending on data format
    if (os.format() == IOstream::ASCII || !contiguous<T>())
    {
        bool uniform = false;

        if (Size > 1 && contiguous<T>())
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
            // Write size (so it is valid dictionary entry) and start delimiter
            os << L.size() << token::BEGIN_BLOCK;

            // Write contents
            os << L[0];

            // Write end delimiter
            os << token::END_BLOCK;
        }
        else if (Size <= 1 ||(Size < 11 && contiguous<T>()))
        {
            // Write start delimiter
            os << token::BEGIN_LIST;

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
            // Write start delimiter
            os << nl << token::BEGIN_LIST;

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
        os.write(reinterpret_cast<const char*>(L.cdata()), Size*sizeof(T));
    }

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, const FixedList&)");

    return os;
}


// ************************************************************************* //
