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

#include "meshToMesh0.H"
#include "volFields.H"
#include "interpolationCellPoint.H"
#include "SubField.H"
#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::meshToMesh0::mapField
(
    Field<Type>& toF,
    const Field<Type>& fromVf,
    const labelList& adr
) const
{
    // Direct mapping of nearest-cell values

    forAll(toF, celli)
    {
        if (adr[celli] != -1)
        {
            toF[celli] = fromVf[adr[celli]];
        }
    }

    // toF.map(fromVf, adr);
}


template<class Type>
void Foam::meshToMesh0::interpolateField
(
    Field<Type>& toF,
    const VolField<Type>& fromVf,
    const labelListList& adr,
    const scalarListList& weights
) const
{
    // Inverse volume weighted interpolation
    forAll(toF, celli)
    {
        const labelList& overlapCells = adr[celli];
        const scalarList& w = weights[celli];

        Type f = Zero;
        forAll(overlapCells, i)
        {
            label fromCelli = overlapCells[i];
            f += fromVf[fromCelli]*w[i];
            toF[celli] = f;
        }
    }
}


template<class Type>
void Foam::meshToMesh0::interpolateField
(
    Field<Type>& toF,
    const VolField<Type>& fromVf,
    const labelList& adr,
    const scalarListList& weights
) const
{
    // Inverse distance weighted interpolation

    // get reference to cellCells
    const labelListList& cc = fromMesh_.cellCells();

    forAll(toF, celli)
    {
        if (adr[celli] != -1)
        {
            const labelList& neighbours = cc[adr[celli]];
            const scalarList& w = weights[celli];

            Type f = fromVf[adr[celli]]*w[0];

            for (label ni = 1; ni < w.size(); ni++)
            {
                f += fromVf[neighbours[ni - 1]]*w[ni];
            }

            toF[celli] = f;
        }
    }
}


template<class Type>
void Foam::meshToMesh0::interpolateField
(
    Field<Type>& toF,
    const VolField<Type>& fromVf,
    const labelList& adr,
    const vectorField& centres
) const
{
    // Cell-Point interpolation
    interpolationCellPoint<Type> interpolator(fromVf);

    forAll(toF, celli)
    {
        if (adr[celli] != -1)
        {
            toF[celli] = interpolator.interpolate(centres[celli], adr[celli]);
        }
    }
}


template<class Type>
void Foam::meshToMesh0::interpolateInternalField
(
    Field<Type>& toF,
    const VolField<Type>& fromVf,
    meshToMesh0::order ord
) const
{
    if (fromVf.mesh() != fromMesh_)
    {
        FatalErrorInFunction
            << "the argument field does not correspond to the right mesh. "
            << "Field size: " << fromVf.size()
            << " mesh size: " << fromMesh_.nCells()
            << exit(FatalError);
    }

    if (toF.size() != toMesh_.nCells())
    {
        FatalErrorInFunction
            << "the argument field does not correspond to the right mesh. "
            << "Field size: " << toF.size()
            << " mesh size: " << toMesh_.nCells()
            << exit(FatalError);
    }

    switch(ord)
    {
        case MAP:
            mapField(toF, fromVf, cellAddressing_);
        break;

        case INTERPOLATE:
        {
            interpolateField
            (
                toF,
                fromVf,
                cellAddressing_,
                inverseDistanceWeights()
            );
            break;
        }
        case CELL_POINT_INTERPOLATE:
        {
            interpolateField
            (
                toF,
                fromVf,
                cellAddressing_,
                toMesh_.cellCentres()
            );

            break;
        }
        case CELL_VOLUME_WEIGHT:
        {
            const labelListList& cellToCell = cellToCellAddressing();
            const scalarListList& invVolWeights = inverseVolumeWeights();

            interpolateField
            (
                toF,
                fromVf,
                cellToCell,
                invVolWeights
            );
            break;
        }
        default:
            FatalErrorInFunction
                << "unknown interpolation scheme " << ord
                << exit(FatalError);
    }
}


template<class Type>
void Foam::meshToMesh0::interpolateInternalField
(
    Field<Type>& toF,
    const tmp<VolField<Type>>& tfromVf,
    meshToMesh0::order ord
) const
{
    interpolateInternalField(toF, tfromVf(), ord);
    tfromVf.clear();
}


template<class Type>
void Foam::meshToMesh0::interpolate
(
    VolField<Type>& toVf,
    const VolField<Type>& fromVf,
    meshToMesh0::order ord
) const
{
    interpolateInternalField(toVf, fromVf, ord);

    typename VolField<Type>::
        Boundary& toVfBf = toVf.boundaryFieldRef();

    forAll(toMesh_.boundaryMesh(), patchi)
    {
        const fvPatch& toPatch = toMesh_.boundary()[patchi];

        if (cuttingPatches_.found(toPatch.name()))
        {
            switch(ord)
            {
                case MAP:
                {
                    mapField
                    (
                        toVfBf[patchi],
                        fromVf,
                        boundaryAddressing_[patchi]
                    );
                    break;
                }

                case INTERPOLATE:
                {
                    interpolateField
                    (
                        toVfBf[patchi],
                        fromVf,
                        boundaryAddressing_[patchi],
                        toPatch.Cf()
                    );
                    break;
                }

                case CELL_POINT_INTERPOLATE:
                {
                    interpolateField
                    (
                        toVfBf[patchi],
                        fromVf,
                        boundaryAddressing_[patchi],
                        toPatch.Cf()
                    );
                    break;
                }
                case CELL_VOLUME_WEIGHT:
                {
                    break;
                }

                default:
                    FatalErrorInFunction
                        << "unknown interpolation scheme " << ord
                        << exit(FatalError);
            }

            if (isA<mixedFvPatchField<Type>>(toVfBf[patchi]))
            {
                refCast<mixedFvPatchField<Type>>
                (
                    toVfBf[patchi]
                ).refValue() = toVfBf[patchi];
            }
        }
        else if
        (
            patchMap_.found(toPatch.name())
         && fromMeshPatches_.found(patchMap_.find(toPatch.name())())
        )
        {
            mapField
            (
                toVfBf[patchi],
                fromVf.boundaryField()
                [
                    fromMeshPatches_.find(patchMap_.find(toPatch.name())())()
                ],
                boundaryAddressing_[patchi]
            );
        }
    }
}


template<class Type>
void Foam::meshToMesh0::interpolate
(
    VolField<Type>& toVf,
    const tmp<VolField<Type>>& tfromVf,
    meshToMesh0::order ord
) const
{
    interpolate(toVf, tfromVf(), ord);
    tfromVf.clear();
}


template<class Type>
Foam::tmp<Foam::VolField<Type>>
Foam::meshToMesh0::interpolate
(
    const VolField<Type>& fromVf,
    meshToMesh0::order ord
) const
{
    // Create and map the internal-field values
    Field<Type> internalField(toMesh_.nCells());
    interpolateInternalField(internalField, fromVf, ord);

    // check whether both meshes have got the same number
    // of boundary patches
    if (fromMesh_.boundary().size() != toMesh_.boundary().size())
    {
        FatalErrorInFunction
            << "Incompatible meshes: different number of boundaries, "
               "only internal field may be interpolated"
            << exit(FatalError);
    }

    // Create and map the patch field values
    PtrList<fvPatchField<Type>> patchFields
    (
        boundaryAddressing_.size()
    );

    forAll(boundaryAddressing_, patchi)
    {
        patchFields.set
        (
            patchi,
            fvPatchField<Type>::New
            (
                fromVf.boundaryField()[patchi],
                toMesh_.boundary()[patchi],
                DimensionedField<Type, volMesh>::null(),
                patchFieldInterpolator
                (
                    boundaryAddressing_[patchi]
                )
            )
        );
    }


    // Create the complete field from the pieces
    tmp<VolField<Type>> ttoF
    (
        new VolField<Type>
        (
            IOobject
            (
                "interpolated(" + fromVf.name() + ')',
                toMesh_.time().name(),
                toMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            toMesh_,
            fromVf.dimensions(),
            internalField,
            patchFields
        )
    );

    return ttoF;
}


template<class Type>
Foam::tmp<Foam::VolField<Type>>
Foam::meshToMesh0::interpolate
(
    const tmp<VolField<Type>>& tfromVf,
    meshToMesh0::order ord
) const
{
    tmp<VolField<Type>> tint =
        interpolate(tfromVf(), ord);
    tfromVf.clear();

    return tint;
}


// ************************************************************************* //
