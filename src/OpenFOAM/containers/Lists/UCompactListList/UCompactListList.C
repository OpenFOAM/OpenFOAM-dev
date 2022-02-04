/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "UCompactListList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
Foam::labelList Foam::UCompactListList<T>::sizes() const
{
    labelList rowSizes(size());

    if (rowSizes.size() > 0)
    {
        forAll(rowSizes, i)
        {
            rowSizes[i] = offsets_[i+1] - offsets_[i];
        }
    }
    return rowSizes;
}


template<class T>
template<class Container>
Foam::List<Container> Foam::UCompactListList<T>::list() const
{
    List<Container> ll(size());

    forAll(ll, i)
    {
        ll[i] = Container(operator[](i));
    }

    return ll;
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "UCompactListListIO.C"

// ************************************************************************* //
