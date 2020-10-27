/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2020 OpenFOAM Foundation
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type, class FOType>
bool Foam::functionObjects::fieldsExpression::calcFieldTypes(FOType& fo)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    if (foundObject<VolFieldType>(fieldNames_[0]))
    {
        return store
        (
            resultName_,
            fo.template calcFieldType<VolFieldType>()
        );
    }
    else if (foundObject<SurfaceFieldType>(fieldNames_[0]))
    {
        return store
        (
            resultName_,
            fo.template calcFieldType<SurfaceFieldType>()
        );
    }
    else
    {
        return false;
    }
}


template<class Type, class FOType>
bool Foam::functionObjects::fieldsExpression::calcType(FOType& fo)
{
    return calcFieldTypes<Type>(fo);
}


template<class FOType>
bool Foam::functionObjects::fieldsExpression::calcAllTypes(FOType& fo)
{
    bool processed = false;

    #define processType(fieldType, none)                                       \
        processed = processed || fo.template calcType<fieldType>(fo);
    FOR_ALL_FIELD_TYPES(processType)

    return processed;
}


// ************************************************************************* //
