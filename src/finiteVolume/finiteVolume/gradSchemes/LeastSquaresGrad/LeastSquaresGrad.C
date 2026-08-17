/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2026 OpenFOAM Foundation
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

#include "LeastSquaresGrad.H"
#include "LeastSquaresVectors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class Stencil>
void Foam::fv::LeastSquaresGrad<Type, Stencil>::calcGrad
(
    VolInternalField<typename outerProduct<vector, Type>::type>& grad,
    const VolField<Type>& vtf
) const
{
    const fvMesh& mesh = vtf.mesh();

    // Get reference to least square vectors
    const LeastSquaresVectors<Stencil>& lsv = LeastSquaresVectors<Stencil>::New
    (
        mesh
    );

    const extendedCentredCellToCellStencil& stencil = lsv.stencil();
    const List<List<label>>& stencilAddr = stencil.stencil();
    const List<List<vector>>& lsvs = lsv.vectors();

    // Construct flat version of vtf
    // including all values referred to by the stencil
    List<Type> flatVtf(stencil.map().constructSize(), Zero);

    // Insert internal values
    forAll(vtf, celli)
    {
        flatVtf[celli] = vtf[celli];
    }

    // Insert boundary values
    forAll(vtf.boundaryField(), patchi)
    {
        const fvPatchField<Type>& ptf = vtf.boundaryField()[patchi];

        label nCompact =
            ptf.patch().start()
          - mesh.nInternalFaces()
          + mesh.nCells();

        forAll(ptf, i)
        {
            flatVtf[nCompact++] = ptf[i];
        }
    }

    // Do all swapping to complete flatVtf
    stencil.map().distribute(flatVtf);

    grad = Zero;

    // Accumulate the cell-centred gradient from the
    // weighted least-squares vectors and the flattened field values
    forAll(stencilAddr, celli)
    {
        const labelList& compactCells = stencilAddr[celli];
        const List<vector>& lsvc = lsvs[celli];

        forAll(compactCells, i)
        {
            grad[celli] += lsvc[i]*flatVtf[compactCells[i]];
        }
    }
}


// ************************************************************************* //
