/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::blockMeshTools::read
(
    Istream& is,
    List<T>& L,
    const dictionary& dict
)
{
    token firstToken(is);

    if (firstToken.isLabel())
    {
        label s = firstToken.labelToken();

        // Set list length to that read
        L.setSize(s);

        // Read list contents depending on data format

        // Read beginning of contents
        char delimiter = is.readBeginList("List");

        if (s)
        {
            if (delimiter == token::BEGIN_LIST)
            {
                for (label i=0; i<s; i++)
                {
                    read(is, L[i], dict);
                }
            }
        }

        // Read end of contents
        is.readEndList("List");
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

        SLList<T> sll;

        while (true)
        {
            token t(is);
            if (t.isPunctuation() && t.pToken() == token::END_LIST)
            {
                break;
            }
            is.putBack(t);
            T elem;
            read(is, elem, dict);
            sll.append(elem);
        }

        // Convert the singly-linked list to this list
        L = sll;
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int> or '(', found "
            << firstToken.info()
            << exit(FatalIOError);
    }
}


template<class T>
Foam::List<T> Foam::blockMeshTools::read
(
    Istream& is,
    const dictionary& dict
)
{
    List<T> L;
    read(is, L, dict);
    return L;
}


// ************************************************************************* //
