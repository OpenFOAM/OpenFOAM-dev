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

\*---------------------------------------------------------------------------*/

#include "extendedCellToFaceStencil.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::extendedUpwindCellToFaceStencil::weightedSum
(
    const surfaceScalarField& phi,
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    const List<List<scalar>>& ownWeights,
    const List<List<scalar>>& neiWeights
) const
{
    const fvMesh& mesh = fld.mesh();

    // Collect internal and boundary values
    List<List<Type>> ownFld;
    collectData(ownMap(), ownStencil(), fld, ownFld);
    List<List<Type>> neiFld;
    collectData(neiMap(), neiStencil(), fld, neiFld);

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
        if (phi[facei] > 0)
        {
            // Flux out of owner. Use upwind (= owner side) stencil.
            const List<Type>& stField = ownFld[facei];
            const List<scalar>& stWeight = ownWeights[facei];

            forAll(stField, i)
            {
                sf[facei] += stField[i]*stWeight[i];
            }
        }
        else
        {
            const List<Type>& stField = neiFld[facei];
            const List<scalar>& stWeight = neiWeights[facei];

            forAll(stField, i)
            {
                sf[facei] += stField[i]*stWeight[i];
            }
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
                if (phi.boundaryField()[patchi][i] > 0)
                {
                    // Flux out of owner. Use upwind (= owner side) stencil.
                    const List<Type>& stField = ownFld[facei];
                    const List<scalar>& stWeight = ownWeights[facei];

                    forAll(stField, j)
                    {
                        pSfCorr[i] += stField[j]*stWeight[j];
                    }
                }
                else
                {
                    const List<Type>& stField = neiFld[facei];
                    const List<scalar>& stWeight = neiWeights[facei];

                    forAll(stField, j)
                    {
                        pSfCorr[i] += stField[j]*stWeight[j];
                    }
                }
                facei++;
            }
        }
    }

    return tsfCorr;
}


// ************************************************************************* //
