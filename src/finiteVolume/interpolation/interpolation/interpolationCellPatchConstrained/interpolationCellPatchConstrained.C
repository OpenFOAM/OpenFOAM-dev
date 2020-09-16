/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "interpolationCellPatchConstrained.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationCellPatchConstrained<Type>::interpolationCellPatchConstrained
(
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
:
    fieldInterpolation<Type, interpolationCellPatchConstrained<Type>>(psi)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::interpolationCellPatchConstrained<Type>::interpolate
(
    const vector& pt,
    const label celli,
    const label facei
) const
{
    if (facei >= 0 && facei >= this->psi_.mesh().nInternalFaces())
    {
        // Use boundary value
        const polyBoundaryMesh& pbm = this->psi_.mesh().boundaryMesh();
        label patchi = pbm.patchID()[facei-this->psi_.mesh().nInternalFaces()];
        label patchFacei = pbm[patchi].whichFace(facei);

        return this->psi_.boundaryField()[patchi][patchFacei];
    }
    else
    {
        return this->psi_[celli];
    }
}


// ************************************************************************* //
