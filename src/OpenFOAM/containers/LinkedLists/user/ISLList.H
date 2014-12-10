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

Class
    Foam::ISLList

Description
    Intrusive singly-linked list.

\*---------------------------------------------------------------------------*/

#ifndef ISLList_H
#define ISLList_H

#include "ILList.H"
#include "SLListBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ISLList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class ISLList
:
    public ILList<SLListBase, T>
{

public:

    // Constructors

        //- Null construct
        ISLList()
        {}

        //- Construct given initial T
        ISLList(const T& a)
        :
            ILList<SLListBase, T>(a)
        {}

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        ISLList(Istream& is, const INew& inewt)
        :
            ILList<SLListBase, T>(is, inewt)
        {}

        //- Construct from Istream
        ISLList(Istream& is)
        :
            ILList<SLListBase, T>(is)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
