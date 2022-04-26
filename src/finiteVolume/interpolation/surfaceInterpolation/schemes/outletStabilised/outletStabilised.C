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

#include "outletStabilised.H"
#include "zeroGradientFvPatchField.H"
#include "mixedFvPatchField.H"
#include "directionMixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline Foam::tmp<Foam::surfaceScalarField>
Foam::outletStabilised<Type>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<surfaceScalarField> tw = tScheme_().weights(vf);
    surfaceScalarField& w = tw.ref();

    const fvMesh& mesh_ = this->mesh();
    const cellList& cells = mesh_.cells();

    forAll(vf.boundaryField(), patchi)
    {
        if
        (
            isA<zeroGradientFvPatchField<Type>>
                (vf.boundaryField()[patchi])
         || isA<mixedFvPatchField<Type>>(vf.boundaryField()[patchi])
         || isA<directionMixedFvPatchField<Type>>
            (vf.boundaryField()[patchi])
        )
        {
            const labelList& pFaceCells =
                mesh_.boundary()[patchi].faceCells();

            forAll(pFaceCells, pFacei)
            {
                const cell& pFaceCell = cells[pFaceCells[pFacei]];

                forAll(pFaceCell, fi)
                {
                    label facei = pFaceCell[fi];

                    if (mesh_.isInternalFace(facei))
                    {
                        // Apply upwind interpolation
                        w[facei] = pos0(faceFlux_[facei]);
                    }
                }
            }
        }
    }

    return tw;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
inline Foam::outletStabilised<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (tScheme_().corrected())
    {
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tcorr =
            tScheme_().correction(vf);

        GeometricField<Type, fvsPatchField, surfaceMesh>& corr =
            tcorr.ref();

        const fvMesh& mesh_ = this->mesh();
        const cellList& cells = mesh_.cells();

        forAll(vf.boundaryField(), patchi)
        {
            if
            (
                isA<zeroGradientFvPatchField<Type>>
                    (vf.boundaryField()[patchi])
             || isA<mixedFvPatchField<Type>>
                    (vf.boundaryField()[patchi])
            )
            {
                const labelList& pFaceCells =
                    mesh_.boundary()[patchi].faceCells();

                forAll(pFaceCells, pFacei)
                {
                    const cell& pFaceCell = cells[pFaceCells[pFacei]];

                    forAll(pFaceCell, fi)
                    {
                        label facei = pFaceCell[fi];

                        if (mesh_.isInternalFace(facei))
                        {
                            // Remove correction
                            corr[facei] = Zero;
                        }
                    }
                }
            }
        }

        return tcorr;
    }
    else
    {
        return tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        (
            nullptr
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceInterpolationScheme(outletStabilised);
}

// ************************************************************************* //
