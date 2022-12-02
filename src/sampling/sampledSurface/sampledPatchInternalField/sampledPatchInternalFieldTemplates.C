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

#include "sampledPatchInternalField.H"
#include "interpolationCellPoint.H"
#include "PrimitivePatchInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledSurfaces::patchInternalField::sampleField
(
    const VolField<Type>& vField
) const
{
    // One value per face
    tmp<Field<Type>> tvalues(new Field<Type>(patchFaceLabels().size()));
    Field<Type>& values = tvalues.ref();

    forAll(patchStart(), i)
    {
        // Get patchface wise data by sampling internal field
        Field<Type> interpVals = vField.primitiveField();
        mappers_[i].map().distribute(interpVals);

        // Store at correct position in values
        label end =
        (
            i < patchStart().size()-1
          ? patchStart()[i+1]
          : patchFaceLabels().size()
        );

        for (label triI = patchStart()[i]; triI < end; triI++)
        {
            values[triI] = interpVals[patchFaceLabels()[triI]];
        }
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledSurfaces::patchInternalField::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    label sz = 0;
    forAll(patchIDs(), i)
    {
        sz += mesh().boundaryMesh()[patchIDs()[i]].size();
    }

    Field<Type> allPatchVals(sz);
    sz = 0;

    forAll(patchIDs(), i)
    {
        // Cells on which samples are generated
        const labelList& sampleCells = mappers_[i].cellIndices();

        // Send the sample points to the cells
        const distributionMap& distMap = mappers_[i].map();
        pointField samplePoints(mappers_[i].samplePoints());
        distMap.reverseDistribute(sampleCells.size(), samplePoints);

        // Interpolate values
        Field<Type> sampleValues(sampleCells.size());
        forAll(sampleCells, i)
        {
            if (sampleCells[i] != -1)
            {
                sampleValues[i] =
                    interpolator.interpolate(samplePoints[i], sampleCells[i]);
            }
        }

        // Send the values back to the patch
        Field<Type> patchValues(sampleValues);
        distMap.distribute(patchValues);

        // Insert into the full value field
        SubList<Type>(allPatchVals, patchValues.size(), sz) = patchValues;
        sz += patchValues.size();
    }

    // Interpolate to points. Reconstruct the patch of all faces to aid
    // interpolation.

    labelList meshFaceLabels(allPatchVals.size());
    sz = 0;
    forAll(patchIDs(), i)
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchIDs()[i]];
        forAll(pp, i)
        {
            meshFaceLabels[sz++] = pp.start()+i;
        }
    }

    indirectPrimitivePatch allPatches
    (
        IndirectList<face>(mesh().faces(), meshFaceLabels),
        mesh().points()
    );

    return PrimitivePatchInterpolation<indirectPrimitivePatch>
    (
        allPatches
    ).faceToPointInterpolate(allPatchVals);
}


// ************************************************************************* //
