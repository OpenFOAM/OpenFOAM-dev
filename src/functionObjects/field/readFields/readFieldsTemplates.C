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

#include "readFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::readFields::loadField
(
    const word& fieldName,
    PtrList<VolField<Type>>& vflds,
    PtrList<SurfaceField<Type>>& sflds
) const
{
    if (obr_.foundObject<VolField<Type>>(fieldName))
    {
        DebugInfo
            << "readFields : Field " << fieldName << " already in database"
            << endl;
    }
    else if (obr_.foundObject<SurfaceField<Type>>(fieldName))
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
            time_.name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if
        (
            fieldHeader.headerOk()
         && fieldHeader.headerClassName() == VolField<Type>::typeName
        )
        {
            // Store field locally
            Log << "    Reading " << fieldName << endl;

            label sz = vflds.size();
            vflds.setSize(sz+1);
            vflds.set(sz, new VolField<Type>(fieldHeader, mesh_));
        }
        else if
        (
            fieldHeader.headerOk()
         && fieldHeader.headerClassName() == SurfaceField<Type>::typeName
        )
        {
            // Store field locally
            Log << "    Reading " << fieldName << endl;

            label sz = sflds.size();
            sflds.setSize(sz+1);
            sflds.set(sz, new SurfaceField<Type>(fieldHeader, mesh_));
        }
    }
}


// ************************************************************************* //
