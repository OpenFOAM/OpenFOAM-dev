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

#include "fvPatchDistWave.H"
#include "nonConformalFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::List<Foam::labelPair> Foam::fvPatchDistWave::getChangedPatchAndFaces
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction
)
{
    label nChangedFaces = 0;
    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        nChangedFaces += mesh.boundary()[iter.key()].size();
    }

    List<labelPair> changedPatchAndFaces(nChangedFaces);

    label changedFacei = 0;
    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        const label patchi = iter.key();
        const fvPatch& patch = mesh.boundary()[patchi];

        const bool conformal = !isA<nonConformalFvPatch>(patch);

        forAll(patch, patchFacei)
        {
            // If this patch is conformal then discard faces which have been
            // covered up by non-conformal faces to within the minFaceFraction
            // tolerance
            if (conformal)
            {
                const label polyFacei = mesh.polyFacesBf()[patchi][patchFacei];

                const scalar faceFraction =
                    patch.magSf()[patchFacei]/mesh.magFaceAreas()[polyFacei];

                if (faceFraction < minFaceFraction) continue;
            }

            changedPatchAndFaces[changedFacei] = labelPair(patchi, patchFacei);
            changedFacei++;
        }
    }

    changedPatchAndFaces.resize(changedFacei);

    return changedPatchAndFaces;
}


// ************************************************************************* //
