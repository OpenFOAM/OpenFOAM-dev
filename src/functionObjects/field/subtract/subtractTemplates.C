/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
Foam::tmp<GeoFieldType>
Foam::functionObjects::subtract::subtractFields() const
{
    tmp<GeoFieldType> tresult
    (
        lookupObject<GeoFieldType>(fieldNames_[0])
      - lookupObject<GeoFieldType>(fieldNames_[1])
    );

    for (label i=2; i<fieldNames_.size(); i++)
    {
        tresult = tresult - lookupObject<GeoFieldType>(fieldNames_[i]);
    }

    return tresult;
}


template<class Type>
bool Foam::functionObjects::subtract::calcSubtract()
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    if (foundObject<VolFieldType>(fieldNames_[0]))
    {
        return store
        (
            resultName_,
            subtractFields<VolFieldType>()
        );
    }
    else if (foundObject<SurfaceFieldType>(fieldNames_[0]))
    {
        return store
        (
            resultName_,
            subtractFields<SurfaceFieldType>()
        );
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
