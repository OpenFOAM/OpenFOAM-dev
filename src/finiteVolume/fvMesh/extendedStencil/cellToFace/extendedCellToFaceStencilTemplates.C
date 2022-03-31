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

#include "extendedCellToFaceStencil.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::extendedCellToFaceStencil::collectData
(
    const distributionMap& map,
    const labelListList& stencil,
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    List<List<Type>>& stencilFld
)
{
    // 1. Construct cell data in compact addressing
    List<Type> flatFld(map.constructSize(), Zero);

    // Insert my internal values
    forAll(fld, celli)
    {
        flatFld[celli] = fld[celli];
    }
    // Insert my boundary values
    forAll(fld.boundaryField(), patchi)
    {
        const fvPatchField<Type>& pfld = fld.boundaryField()[patchi];

        label nCompact =
            pfld.patch().start()
           -fld.mesh().nInternalFaces()
           +fld.mesh().nCells();

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
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::extendedCellToFaceStencil::weightedSum
(
    const distributionMap& map,
    const labelListList& stencil,
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    const List<List<scalar>>& stencilWeights
)
{
    const fvMesh& mesh = fld.mesh();

    // Collect internal and boundary values
    List<List<Type>> stencilFld;
    collectData(map, stencil, fld, stencilFld);

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfCorr
    (
        GeometricField<Type, fvsPatchField, surfaceMesh>::New
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
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsfCorr.ref();

    // Internal faces
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        const List<Type>& stField = stencilFld[facei];
        const List<scalar>& stWeight = stencilWeights[facei];

        forAll(stField, i)
        {
            sf[facei] += stField[i]*stWeight[i];
        }
    }

    // Boundaries. Either constrained or calculated so assign value
    // directly (instead of nicely using operator==)
    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
        Boundary& bSfCorr = sf.boundaryFieldRef();

    forAll(bSfCorr, patchi)
    {
        fvsPatchField<Type>& pSfCorr = bSfCorr[patchi];

        if (pSfCorr.coupled())
        {
            label facei = pSfCorr.patch().start();

            forAll(pSfCorr, i)
            {
                const List<Type>& stField = stencilFld[facei];
                const List<scalar>& stWeight = stencilWeights[facei];

                forAll(stField, j)
                {
                    pSfCorr[i] += stField[j]*stWeight[j];
                }

                facei++;
            }
        }
    }

    return tsfCorr;
}


// ************************************************************************* //
