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

#include "nearWallFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::nearWallFields::createFields
(
    PtrList<VolField<Type>>& sfields
) const
{
    const UPtrList<VolField<Type>> fields(mesh_.fields<VolField<Type>>());

    forAll(fields, i)
    {
        const VolField<Type>& field = fields[i];

        if (fieldMap_.found(field.name()))
        {
            const word& sampleFieldName = fieldMap_[field.name()];

            if (obr_.found(sampleFieldName))
            {
                Log << "    a field " << sampleFieldName
                    << " already exists on the mesh."
                    << endl;
            }
            else
            {
                label sz = sfields.size();
                sfields.setSize(sz+1);

                sfields.set
                (
                    sz,
                    new VolField<Type>
                    (
                        IOobject
                        (
                            sampleFieldName,
                            time_.name(),
                            mesh_
                        ),
                        field,
                        calculatedFvPatchScalarField::typeName
                    )
                );

                Log << "    created " << sfields[sz].name()
                    << " to sample " << field.name() << endl;
            }
        }
    }
}


template<class Type>
void Foam::functionObjects::nearWallFields::sampleBoundaryField
(
    const interpolationCellPoint<Type>& interpolator,
    VolField<Type>& field
) const
{
    // Construct flat fields for all patch faces to be sampled
    Field<Type> sampledValues(getPatchDataMapPtr_().constructSize());

    forAll(cellToWalls_, celli)
    {
        const labelList& cData = cellToWalls_[celli];

        forAll(cData, i)
        {
            const point& samplePt = cellToSamples_[celli][i];
            sampledValues[cData[i]] = interpolator.interpolate(samplePt, celli);
        }
    }

    // Send back sampled values to patch faces
    getPatchDataMapPtr_().reverseDistribute
    (
        getPatchDataMapPtr_().constructSize(),
        sampledValues
    );

    typename VolField<Type>::
        Boundary& fieldBf = field.boundaryFieldRef();

    // Pick up data
    label nPatchFaces = 0;
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();

        fvPatchField<Type>& pfield = fieldBf[patchi];

        Field<Type> newField(pfield.size());
        forAll(pfield, i)
        {
            newField[i] = sampledValues[nPatchFaces++];
        }

        pfield == newField;
    }
}


template<class Type>
void Foam::functionObjects::nearWallFields::sampleFields
(
    PtrList<VolField<Type>>& sfields
) const
{
    forAll(sfields, i)
    {
        const word& fieldName = reverseFieldMap_[sfields[i].name()];
        const VolField<Type>& field =
            obr_.lookupObject<VolField<Type>>(fieldName);

        // Take over internal and boundary values
        sfields[i] == field;

        // Construct interpolation method
        interpolationCellPoint<Type> interpolator(field);

        // Override sampled values
        sampleBoundaryField(interpolator, sfields[i]);
    }
}


// ************************************************************************* //
