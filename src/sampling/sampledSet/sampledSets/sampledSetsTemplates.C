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

#include "sampledSets.H"
#include "volFields.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::PtrList<Foam::Field<Type>>
Foam::functionObjects::sampledSets::sampleLocalType
(
    const label seti,
    const wordList& fieldNames,
    HashPtrTable<interpolation<Type>>& interpolations
)
{
    PtrList<Field<Type>> fieldTypeValues(fieldNames.size());

    const sampledSet& s = operator[](seti);

    forAll(fieldNames, fieldi)
    {
        const word& name = fieldNames[fieldi];

        if (mesh_.foundObject<VolField<Type>>(name))
        {
            const VolField<Type>& vf =
                mesh_.lookupObject<VolField<Type>>(name);

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

            const interpolation<Type>& interp = *interpolations[name];

            tmp<Field<Type>> tfield(new Field<Type>(s.size()));
            Field<Type>& field = tfield.ref();

            forAll(s, i)
            {
                const point& position = s.positions()[i];
                const label celli = s.cells()[i];
                const label facei = s.faces()[i];

                if (celli == -1 && facei == -1)
                {
                    // Special condition for illegal sampling points
                    field[i] = pTraits<Type>::max;
                }
                else
                {
                    field[i] = interp.interpolate(position, celli, facei);
                }
            }

            fieldTypeValues.set(fieldi, tfield.ptr());
        }
    }

    return fieldTypeValues;
}


template<class Type>
Foam::PtrList<Foam::Field<Type>>
Foam::functionObjects::sampledSets::sampleType
(
    const label seti,
    const wordList& fieldNames,
    HashPtrTable<interpolation<Type>>& interpolations
)
{
    PtrList<Field<Type>> fieldTypeValues =
        sampleLocalType<Type>(seti, fieldNames, interpolations);

    if (Pstream::parRun())
    {
        forAll(fieldNames, fieldi)
        {
            if (fieldTypeValues.set(fieldi))
            {
                fieldTypeValues[fieldi] =
                    coordSet::gather
                    (
                        fieldTypeValues[fieldi],
                        masterSetOrders_[seti]
                    );
            }
        }
    }

    return fieldTypeValues;
}


// ************************************************************************* //
