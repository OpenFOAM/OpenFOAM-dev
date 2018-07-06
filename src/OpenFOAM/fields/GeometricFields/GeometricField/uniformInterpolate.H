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

#include "GeometricField.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

//- Interpolate selected fields (given by indices and corresponding weights)
//  (boundary type becomes calculated). Fields stored per index. Field gets name
//  "uniformInterpolate(" + fld.name() + ')'.
template<class GeoField>
tmp<GeoField> uniformInterpolate
(
    const HashPtrTable<GeoField, label, Hash<label>>& fields,
    const labelList& indices,
    const scalarField& weights
);

//- Interpolate fields. fieldsCache contains per timeName all loaded fields.
//  Resulting field gets properties according to fieldIO
template<class GeoField>
tmp<GeoField> uniformInterpolate
(
    const IOobject& fieldIO,
    const word& fieldName,
    const wordList& times,
    const scalarField& weights,
    const objectRegistry& fieldsCache
);

//- Interpolate fields. fieldsCache contains per timeName all loaded fields.
//  Resulting field gets properties according to fieldIO
template<class GeoField>
tmp<GeoField> uniformInterpolate
(
    const IOobject& fieldIO,
    const word& fieldName,
    const wordList& times,
    const scalarField& weights,
    const word& registryName = "fieldsCache"
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "uniformInterpolate.C"
#endif

// ************************************************************************* //
