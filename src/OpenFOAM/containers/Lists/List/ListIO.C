/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "List.H"
#include "Istream.H"
#include "token.H"
#include "SLList.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// Construct from Istream
template<class T>
Foam::List<T>::List(Istream& is)
:
    UList<T>(NULL, 0)
{
    operator>>(is, *this);
}


template<class T>
Foam::Istream& Foam::operator>>(Istream& is, List<T>& L)
{
    // Anull list
    L.setSize(0);

    is.fatalCheck("operator>>(Istream&, List<T>&)");

    token firstToken(is);

    is.fatalCheck("operator>>(Istream&, List<T>&) : reading first token");

    if (firstToken.isCompound())
    {
        L.transfer
        (
            dynamicCast<token::Compound<List<T> > >
            (
                firstToken.transferCompoundToken(is)
            )
        );
    }
    else if (firstToken.isLabel())
    {
        label s = firstToken.labelToken();

        // Set list length to that read
        L.setSize(s);

        // Read list contents depending on data format

        if (is.format() == IOstream::ASCII || !contiguous<T>())
        {
            // Read beginning of contents
            char delimiter = is.readBeginList("List");

            if (s)
            {
                if (delimiter == token::BEGIN_LIST)
                {
                    for (register label i=0; i<s; i++)
                    {
                        is >> L[i];

                        is.fatalCheck
                        (
                            "operator>>(Istream&, List<T>&) : reading entry"
                        );
                    }
                }
                else
                {
                    T element;
                    is >> element;

                    is.fatalCheck
                    (
                        "operator>>(Istream&, List<T>&) : "
                        "reading the single entry"
                    );

                    for (register label i=0; i<s; i++)
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
                    "operator>>(Istream&, List<T>&) : reading the binary block"
                );
            }
        }
    }
    else if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorIn("operator>>(Istream&, List<T>&)", is)
                << "incorrect first token, expected '(', found "
                << firstToken.info()
                << exit(FatalIOError);
        }

        // Putback the opening bracket
        is.putBack(firstToken);

        // Now read as a singly-linked list
        SLList<T> sll(is);

        // Convert the singly-linked list to this list
        L = sll;
    }
    else
    {
        FatalIOErrorIn("operator>>(Istream&, List<T>&)", is)
            << "incorrect first token, expected <int> or '(', found "
            << firstToken.info()
            << exit(FatalIOError);
    }

    return is;
}


template<class T>
Foam::List<T> Foam::readList(Istream& is)
{
    List<T> L;
    token firstToken(is);
    is.putBack(firstToken);

    if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorIn("readList<T>(Istream&)", is)
                << "incorrect first token, expected '(', found "
                << firstToken.info()
                << exit(FatalIOError);
        }

        // read via a singly-linked list
        L = SLList<T>(is);
    }
    else
    {
        // create list with a single item
        L.setSize(1);

        is >> L[0];
    }

    return L;
}


// ************************************************************************* //
