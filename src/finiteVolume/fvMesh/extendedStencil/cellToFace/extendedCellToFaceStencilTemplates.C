/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    const mapDistribute& map,
    const labelListList& stencil,
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    List<List<Type> >& stencilFld
)
{
    // 1. Construct cell data in compact addressing
    List<Type> flatFld(map.constructSize(), pTraits<Type>::zero);

    // Insert my internal values
    forAll(fld, cellI)
    {
        flatFld[cellI] = fld[cellI];
    }
    // Insert my boundary values
    forAll(fld.boundaryField(), patchI)
    {
        const fvPatchField<Type>& pfld = fld.boundaryField()[patchI];

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

    forAll(stencil, faceI)
    {
        const labelList& compactCells = stencil[faceI];

        stencilFld[faceI].setSize(compactCells.size());

        forAll(compactCells, i)
        {
            stencilFld[faceI][i] = flatFld[compactCells[i]];
        }
    }
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::extendedCellToFaceStencil::weightedSum
(
    const mapDistribute& map,
    const labelListList& stencil,
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    const List<List<scalar> >& stencilWeights
)
{
    const fvMesh& mesh = fld.mesh();

    // Collect internal and boundary values
    List<List<Type> > stencilFld;
    collectData(map, stencil, fld, stencilFld);

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                fld.name(),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensioned<Type>
            (
                fld.name(),
                fld.dimensions(),
                pTraits<Type>::zero
            )
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsfCorr();

    // Internal faces
    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        const List<Type>& stField = stencilFld[faceI];
        const List<scalar>& stWeight = stencilWeights[faceI];

        forAll(stField, i)
        {
            sf[faceI] += stField[i]*stWeight[i];
        }
    }

    // Boundaries. Either constrained or calculated so assign value
    // directly (instead of nicely using operator==)
    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
        GeometricBoundaryField& bSfCorr = sf.boundaryField();

    forAll(bSfCorr, patchi)
    {
        fvsPatchField<Type>& pSfCorr = bSfCorr[patchi];

        if (pSfCorr.coupled())
        {
            label faceI = pSfCorr.patch().start();

            forAll(pSfCorr, i)
            {
                const List<Type>& stField = stencilFld[faceI];
                const List<scalar>& stWeight = stencilWeights[faceI];

                forAll(stField, j)
                {
                    pSfCorr[i] += stField[j]*stWeight[j];
                }

                faceI++;
            }
        }
    }

    return tsfCorr;
}


// ************************************************************************* //
