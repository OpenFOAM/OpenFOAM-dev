/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "sampledSurfaces.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ListListOps.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::PtrList<Foam::Field<Type>>
Foam::functionObjects::sampledSurfaces::sampleLocalType
(
    const label surfi,
    const wordList& fieldNames,
    HashPtrTable<interpolation<Type>>& interpolations
)
{
    PtrList<Field<Type>> fieldTypeValues(fieldNames.size());

    const sampledSurface& s = operator[](surfi);

    forAll(fieldNames, fieldi)
    {
        const word& name = fieldNames[fieldi];

        if (mesh_.foundObject<VolField<Type>>(name))
        {
            const VolField<Type>& vf =
                mesh_.lookupObject<VolField<Type>>(name);

            if (s.interpolate())
            {
                if (!interpolations.found(name))
                {
                    interpolations.insert
                    (
                        name,
                        interpolation<Type>::New
                        (
                            interpolationScheme_,
                            vf
                        ).ptr()
                    );
                }

                fieldTypeValues.set
                (
                    fieldi,
                    s.interpolate(*interpolations[name]).ptr()
                );
            }
            else
            {
                fieldTypeValues.set(fieldi, s.sample(vf).ptr());
            }
        }
        else if (mesh_.foundObject<SurfaceField<Type>>(name))
        {
            const SurfaceField<Type>& sf =
                mesh_.lookupObject<SurfaceField<Type>>(name);

            fieldTypeValues.set(fieldi, s.sample(sf).ptr());
        }
    }

    return fieldTypeValues;
}


template<class Type>
Foam::PtrList<Foam::Field<Type>>
Foam::functionObjects::sampledSurfaces::sampleType
(
    const label surfi,
    const wordList& fieldNames,
    HashPtrTable<interpolation<Type>>& interpolations
)
{
    // Generate local samples
    PtrList<Field<Type>> fieldTypeValues =
        sampleLocalType<Type>(surfi, fieldNames, interpolations);

    if (Pstream::parRun())
    {
        // Collect values from all processors
        PtrList<List<Field<Type>>> gatheredTypeValues(fieldNames.size());
        forAll(fieldNames, fieldi)
        {
            if (fieldTypeValues.set(fieldi))
            {
                gatheredTypeValues.set
                (
                    fieldi,
                    new List<Field<Type>>(Pstream::nProcs())
                );
                gatheredTypeValues[fieldi][Pstream::myProcNo()] =
                    fieldTypeValues[fieldi];
                Pstream::gatherList(gatheredTypeValues[fieldi]);
            }
        }

        // Clear the local field values
        fieldTypeValues.clear();
        fieldTypeValues.resize(fieldNames.size());

        // Combine on the master
        if (Pstream::master())
        {
            // Combine values into single field
            forAll(fieldNames, fieldi)
            {
                if (gatheredTypeValues.set(fieldi))
                {
                    fieldTypeValues.set
                    (
                        fieldi,
                        new Field<Type>
                        (
                            ListListOps::combine<Field<Type>>
                            (
                                gatheredTypeValues[fieldi],
                                accessOp<Field<Type>>()
                            )
                        )
                    );
                }
            }

            // Renumber point data to correspond to merged points
            forAll(fieldNames, fieldi)
            {
                if (fieldTypeValues.set(fieldi))
                {
                    if
                    (
                        mergeList_[surfi].pointsMap.size()
                     == fieldTypeValues[fieldi].size()
                    )
                    {
                        Field<Type> f(fieldTypeValues[fieldi]);

                        inplaceReorder(mergeList_[surfi].pointsMap, f);
                        f.setSize(mergeList_[surfi].points.size());

                        fieldTypeValues.set(fieldi, new Field<Type>(f, true));
                    }
                }
            }
        }
    }

    return fieldTypeValues;
}


// ************************************************************************* //
