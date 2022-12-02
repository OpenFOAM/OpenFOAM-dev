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

#include "extendedFaceToCellStencil.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::extendedFaceToCellStencil::collectData
(
    const distributionMap& map,
    const labelListList& stencil,
    const SurfaceField<Type>& fld,
    List<List<Type>>& stencilFld
)
{
    // 1. Construct face data in compact addressing
    List<Type> flatFld(map.constructSize(), Zero);

    // Insert my internal values
    forAll(fld, celli)
    {
        flatFld[celli] = fld[celli];
    }
    // Insert my boundary values
    forAll(fld.boundaryField(), patchi)
    {
        const fvsPatchField<Type>& pfld = fld.boundaryField()[patchi];

        label nCompact = pfld.patch().start();

        forAll(pfld, i)
        {
            flatFld[nCompact++] = pfld[i];
        }
    }

    // Do all swapping
    map.distribute(flatFld);

    // 2. Pull to stencil
    stencilFld.setSize(stencil.size());

    forAll(stencil, facei)
    {
        const labelList& compactCells = stencil[facei];

        stencilFld[facei].setSize(compactCells.size());

        forAll(compactCells, i)
        {
            stencilFld[facei][i] = flatFld[compactCells[i]];
        }
    }
}


template<class Type>
Foam::tmp<Foam::VolField<Type>>
Foam::extendedFaceToCellStencil::weightedSum
(
    const distributionMap& map,
    const labelListList& stencil,
    const SurfaceField<Type>& fld,
    const List<List<scalar>>& stencilWeights
)
{
    const fvMesh& mesh = fld.mesh();

    // Collect internal and boundary values
    List<List<Type>> stencilFld;
    collectData(map, stencil, fld, stencilFld);

    tmp<VolField<Type>> tsfCorr
    (
        VolField<Type>::New
        (
            fld.name(),
            mesh,
            dimensioned<Type>
            (
                fld.name(),
                fld.dimensions(),
                Zero
            )
        )
    );
    VolField<Type>& sf = tsfCorr.ref();

    // cells
    forAll(sf, celli)
    {
        const List<Type>& stField = stencilFld[celli];
        const List<scalar>& stWeight = stencilWeights[celli];

        forAll(stField, i)
        {
            sf[celli] += stField[i]*stWeight[i];
        }
    }

    // Boundaries values?

    return tsfCorr;
}


// ************************************************************************* //
