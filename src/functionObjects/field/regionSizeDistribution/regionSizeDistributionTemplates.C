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

#include "regionSizeDistribution.H"
#include "regionSplit.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::Map<Type> Foam::functionObjects::regionSizeDistribution::regionSum
(
    const regionSplit& regions,
    const Field<Type>& fld
) const
{
    // Per region the sum of fld
    Map<Type> regionToSum(regions.nRegions()/Pstream::nProcs());

    forAll(fld, celli)
    {
        label regionI = regions[celli];

        typename Map<Type>::iterator fnd = regionToSum.find(regionI);
        if (fnd == regionToSum.end())
        {
            regionToSum.insert(regionI, fld[celli]);
        }
        else
        {
            fnd() += fld[celli];
        }
    }
    Pstream::mapCombineGather(regionToSum, plusEqOp<Type>());
    Pstream::mapCombineScatter(regionToSum);

    return regionToSum;
}


template<class Type>
Foam::List<Type> Foam::functionObjects::regionSizeDistribution::extractData
(
    const UList<label>& keys,
    const Map<Type>& regionData
) const
{
    List<Type> sortedData(keys.size());

    forAll(keys, i)
    {
        sortedData[i] = regionData[keys[i]];
    }
    return sortedData;
}


// ************************************************************************* //
