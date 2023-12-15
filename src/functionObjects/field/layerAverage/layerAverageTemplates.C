/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "layerAverage.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
T Foam::functionObjects::layerAverage::symmetricCoeff() const
{
    return pTraits<T>::one;
}


template<class T>
Foam::tmp<Foam::Field<T>> Foam::functionObjects::layerAverage::sum
(
    const VolInternalField<T>& cellField
) const
{
    tmp<Field<T>> tlayerField(new Field<T>(nLayers_, Zero));
    Field<T>& layerField = tlayerField.ref();

    forAll(cellLayer_, celli)
    {
        if (cellLayer_[celli] != -1)
        {
            layerField[cellLayer_[celli]] += cellField[celli];
        }
    }

    Pstream::listCombineGather(layerField, plusEqOp<T>());
    Pstream::listCombineScatter(layerField);

    return tlayerField;
}


template<class T>
Foam::tmp<Foam::Field<T>> Foam::functionObjects::layerAverage::average
(
    const tmp<VolInternalField<scalar>>& cellWeight,
    const tmp<Field<scalar>>& layerWeight,
    const VolInternalField<T>& cellField
) const
{
    tmp<Field<T>> tlayerField
    (
        cellWeight.valid()
      ? sum<T>(mesh_.V()*cellWeight*cellField)/layerWeight
      : sum<T>(mesh_.V()*cellField)/layerVolume_
    );

    // Handle symmetry
    if (symmetric_)
    {
        Field<T>& layerField = tlayerField.ref();

        const T coeff = symmetricCoeff<T>();

        for (label i=0; i<nLayers_/2; i++)
        {
            const label j = nLayers_ - i - 1;

            layerField[i] =
                (layerField[i] + cmptMultiply(coeff, layerField[j]))/2;
        }

        layerField.setSize(nLayers_/2);
    }

    return tlayerField;
}


// ************************************************************************* //
