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

Description

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::cell::cell()
{}


inline Foam::cell::cell(label s)
:
    labelList(s, -1)
{}


inline Foam::cell::cell(const labelUList& lst)
:
    labelList(lst)
{}


inline Foam::cell::cell(labelList&& lst)
:
    labelList(move(lst))
{}


inline Foam::cell::cell(Istream& is)
:
    labelList(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::cell::nFaces() const
{
    return size();
}


inline bool Foam::operator!=(const cell& a, const cell& b)
{
    return (!(a == b));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::cell::operator=(labelList&& l)
{
    labelList::operator=(move(l));
}


// ************************************************************************* //
