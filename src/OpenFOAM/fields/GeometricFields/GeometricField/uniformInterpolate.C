/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GeoField>
Foam::tmp<GeoField> Foam::uniformInterpolate
(
    const HashPtrTable<GeoField, label, Hash<label>>& fields,
    const labelList& indices,
    const scalarField& weights
)
{
    const GeoField& field0 = *(*fields.begin());

    // Interpolate
    tmp<GeoField> tfld
    (
        new GeoField
        (
            IOobject
            (
                "uniformInterpolate(" + field0.name() + ')',
                field0.time().timeName(),
                field0.db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            weights[0]*(*fields[indices[0]])
        )
    );
    GeoField& fld = tfld();

    for (label i = 1; i < indices.size(); ++i)
    {
        fld += weights[i]*(*fields[indices[i]]);
    }

    return tfld;
}


template<class GeoField>
Foam::tmp<GeoField> Foam::uniformInterpolate
(
    const IOobject& fieldIO,
    const word& fieldName,
    const wordList& times,
    const scalarField& weights,
    const objectRegistry& fieldsCache
)
{
    // Look up the first field
    const objectRegistry& time0Fields = fieldsCache.lookupObject
    <
        const objectRegistry
    >
    (
        times[0]
    );
    const GeoField& field0 = time0Fields.lookupObject
    <
        const GeoField
    >
    (
        fieldName
    );


    // Interpolate
    tmp<GeoField> tfld(new GeoField(fieldIO, weights[0]*field0));
    GeoField& fld = tfld.ref();

    for (label i = 1; i < times.size(); ++i)
    {
        const objectRegistry& timeIFields = fieldsCache.lookupObject
        <
            const objectRegistry
        >
        (
            times[i]
        );
        const GeoField& fieldi = timeIFields.lookupObject
        <
            const GeoField
        >
        (
            fieldName
        );

        fld += weights[i]*fieldi;
    }

    return tfld;
}


template<class GeoField>
Foam::tmp<GeoField> Foam::uniformInterpolate
(
    const IOobject& fieldIO,
    const word& fieldName,
    const wordList& times,
    const scalarField& weights,
    const word& registryName
)
{
    return uniformInterpolate<GeoField>
    (
        fieldIO,
        fieldName,
        times,
        weights,
        fieldIO.db().subRegistry(registryName, true)
    );
}


// ************************************************************************* //
