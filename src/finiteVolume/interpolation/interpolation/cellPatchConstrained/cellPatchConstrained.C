/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "cellPatchConstrained.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolations::cellPatchConstrained<Type>::cellPatchConstrained
(
    const VolField<Type>& psi
)
:
    fieldInterpolation<Type, cellPatchConstrained<Type>>(psi)
{}


template<class Type>
Foam::interpolations::cellPatchConstrained<Type>::cellPatchConstrained
(
    const cellPatchConstrained<Type>& i
)
:
    fieldInterpolation<Type, cellPatchConstrained<Type>>(i)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::interpolations::cellPatchConstrained<Type>::interpolate
(
    const vector& pt,
    const label celli,
    const label facei
) const
{
    const fvMesh& mesh = this->psi_.mesh();

    if (facei >= 0 && facei >= this->psi_.mesh().nInternalFaces())
    {
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        const label patchi = pbm.patchIndices()[facei - mesh.nInternalFaces()];

        if (!pbm[patchi].empty())
        {
            const label patchFacei = pbm[patchi].whichFace(facei);

            return this->psi_.boundaryField()[patchi][patchFacei];
        }
    }

    return this->psi_[celli];
}


// ************************************************************************* //
