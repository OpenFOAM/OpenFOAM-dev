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

#include "LList.H"
#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::LList<LListBase, T>::LList(Istream& is)
{
    operator>>(is, *this);
}


// * * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::Istream& Foam::operator>>(Istream& is, LList<LListBase, T>& L)
{
    // Anull list
    L.clear();

    is.fatalCheck(" operator>>(Istream&, LList<LListBase, T>&)");

    token firstToken(is);

    is.fatalCheck
    (
        " operator>>(Istream&, LList<LListBase, T>&) : reading first token"
    );

    if (firstToken.isLabel())
    {
        label s = firstToken.labelToken();

        // Read beginning of contents
        char delimiter = is.readBeginList("LList<LListBase, T>");

        if (s)
        {
            if (delimiter == token::BEGIN_LIST)
            {
                for (register label i=0; i<s; ++i)
                {
                    T element;
                    is >> element;
                    L.append(element);
                }
            }
            else
            {
                T element;
                is >> element;

                for (register label i=0; i<s; ++i)
                {
                    L.append(element);
                }
            }
        }

        // Read end of contents
        is.readEndList("LList");
    }
    else if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorIn
            (
                " operator>>(Istream&, LList<LListBase, T>&)",
                is
            )   << "incorrect first token, '(', found " << firstToken.info()
                << exit(FatalIOError);
        }

        token lastToken(is);
        is.fatalCheck(" operator>>(Istream&, LList<LListBase, T>&)");

        while
        (
           !(
                lastToken.isPunctuation()
             && lastToken.pToken() == token::END_LIST
            )
        )
        {
            is.putBack(lastToken);
            T element;
            is >> element;
            L.append(element);

            is >> lastToken;
            is.fatalCheck(" operator>>(Istream&, LList<LListBase, T>&)");
        }
    }
    else
    {
        FatalIOErrorIn(" operator>>(Istream&, LList<LListBase, T>&)", is)
            << "incorrect first token, expected <int> or '(', found "
            << firstToken.info()
            << exit(FatalIOError);
    }

    // Check state of IOstream
    is.fatalCheck(" operator>>(Istream&, LList<LListBase,>&)");

    return is;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::Ostream& Foam::operator<<(Ostream& os, const LList<LListBase, T>& lst)
{
    // Write size
    os << nl << lst.size();

    // Write beginning of contents
    os << nl << token::BEGIN_LIST << nl;

    // Write contents
    for
    (
        typename LList<LListBase, T>::const_iterator iter = lst.begin();
        iter != lst.end();
        ++iter
    )
    {
        os << iter() << nl;
    }

    // Write end of contents
    os << token::END_LIST;

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, const LList<LListBase, T>&)");

    return os;
}


// ************************************************************************* //
