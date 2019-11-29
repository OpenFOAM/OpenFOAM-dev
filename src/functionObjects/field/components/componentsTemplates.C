/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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

#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GeoFieldType>
bool Foam::functionObjects::components::calcFieldComponents()
{
    typedef typename GeoFieldType::value_type Type;

    const GeoFieldType& field(lookupObject<GeoFieldType>(fieldName_));

    resultNames_.setSize(Type::nComponents);

    bool stored = true;

    for (direction i=0; i<Type::nComponents; i++)
    {
        resultNames_[i] = fieldName_ + word(Type::componentNames[i]);
        stored = stored && store(resultNames_[i], field.component(i));
    }

    return stored;
}


template<class Type>
bool Foam::functionObjects::components::calcComponents()
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    if (foundObject<VolFieldType>(fieldName_))
    {
        return calcFieldComponents<VolFieldType>();
    }
    else if (foundObject<SurfaceFieldType>(fieldName_))
    {
        return calcFieldComponents<SurfaceFieldType>();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
