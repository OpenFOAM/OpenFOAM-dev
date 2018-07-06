/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "interpolationWeights.H"
#include "ListOps.H"
#include "IOobject.H"
#include "HashSet.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ListType1, class ListType2>
typename Foam::outerProduct
<
    typename ListType1::value_type,
    typename ListType2::value_type
>::type
Foam::interpolationWeights::weightedSum
(
    const ListType1& f1,
    const ListType2& f2
)
{
    typedef typename outerProduct
    <
        typename ListType1::value_type,
        typename ListType2::value_type
    >::type returnType;

    if (f1.size())
    {
        returnType SumProd = f1[0]*f2[0];
        for (label i = 1; i < f1.size(); ++i)
        {
            SumProd += f1[i]*f2[i];
        }
        return SumProd;
    }
    else
    {
        return Zero;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
