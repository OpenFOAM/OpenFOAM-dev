/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "readFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::readFields::loadField
(
    const word& fieldName,
    PtrList<GeometricField<Type, fvPatchField, volMesh>>& vflds,
    PtrList<GeometricField<Type, fvsPatchField, surfaceMesh>>& sflds
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    if (obr_.foundObject<VolFieldType>(fieldName))
    {
        DebugInfo
            << "readFields : Field " << fieldName << " already in database"
            << endl;
    }
    else if (obr_.foundObject<SurfaceFieldType>(fieldName))
    {
        DebugInfo
            << "readFields : Field " << fieldName
            << " already in database" << endl;
    }
    else
    {
        IOobject fieldHeader
        (
            fieldName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if
        (
            fieldHeader.typeHeaderOk<VolFieldType>(false)
         && fieldHeader.headerClassName() == VolFieldType::typeName
        )
        {
            // Store field locally
            Log << "    Reading " << fieldName << endl;

            label sz = vflds.size();
            vflds.setSize(sz+1);
            vflds.set(sz, new VolFieldType(fieldHeader, mesh_));
        }
        else if
        (
            fieldHeader.typeHeaderOk<SurfaceFieldType>(false)
         && fieldHeader.headerClassName() == SurfaceFieldType::typeName
        )
        {
            // Store field locally
            Log << "    Reading " << fieldName << endl;

            label sz = sflds.size();
            sflds.setSize(sz+1);
            sflds.set(sz, new SurfaceFieldType(fieldHeader, mesh_));
        }
    }
}


// ************************************************************************* //
