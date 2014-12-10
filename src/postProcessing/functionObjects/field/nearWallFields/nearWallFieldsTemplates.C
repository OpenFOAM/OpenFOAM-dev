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

#include "nearWallFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::nearWallFields::createFields
(
    PtrList<GeometricField<Type, fvPatchField, volMesh> >& sflds
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> vfType;

    HashTable<const vfType*> flds(obr_.lookupClass<vfType>());

    forAllConstIter(typename HashTable<const vfType*>, flds, iter)
    {
        const vfType& fld = *iter();

        if (fieldMap_.found(fld.name()))
        {
            const word& sampleFldName = fieldMap_[fld.name()];

            if (obr_.found(sampleFldName))
            {
                Info<< "    a field " << sampleFldName
                    << " already exists on the mesh."
                    << endl;
            }
            else
            {
                label sz = sflds.size();
                sflds.setSize(sz+1);

                IOobject io(fld);
                io.readOpt() = IOobject::NO_READ;
                io.writeOpt() = IOobject::NO_WRITE;
                io.rename(sampleFldName);

                sflds.set(sz, new vfType(io, fld));

                Info<< "    created " << sflds[sz].name() << " to sample "
                    << fld.name() << endl;
            }
        }
    }
}


template<class Type>
void Foam::nearWallFields::sampleBoundaryField
(
    const interpolationCellPoint<Type>& interpolator,
    GeometricField<Type, fvPatchField, volMesh>& fld
) const
{
    // Construct flat fields for all patch faces to be sampled
    Field<Type> sampledValues(getPatchDataMapPtr_().constructSize());

    forAll(cellToWalls_, cellI)
    {
        const labelList& cData = cellToWalls_[cellI];

        forAll(cData, i)
        {
            const point& samplePt = cellToSamples_[cellI][i];
            sampledValues[cData[i]] = interpolator.interpolate(samplePt, cellI);
        }
    }

    // Send back sampled values to patch faces
    getPatchDataMapPtr_().reverseDistribute
    (
        getPatchDataMapPtr_().constructSize(),
        sampledValues
    );

    // Pick up data
    label nPatchFaces = 0;
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchI = iter.key();

        fvPatchField<Type>& pfld = fld.boundaryField()[patchI];

        Field<Type> newFld(pfld.size());
        forAll(pfld, i)
        {
            newFld[i] = sampledValues[nPatchFaces++];
        }

        pfld == newFld;
    }
}


template<class Type>
void Foam::nearWallFields::sampleFields
(
    PtrList<GeometricField<Type, fvPatchField, volMesh> >& sflds
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> vfType;

    forAll(sflds, i)
    {
        const word& fldName = reverseFieldMap_[sflds[i].name()];
        const vfType& fld = obr_.lookupObject<vfType>(fldName);

        // Take over internal and boundary values
        sflds[i] == fld;

        // Construct interpolation method
        interpolationCellPoint<Type> interpolator(fld);

        // Override sampled values
        sampleBoundaryField(interpolator, sflds[i]);
    }
}


// ************************************************************************* //
