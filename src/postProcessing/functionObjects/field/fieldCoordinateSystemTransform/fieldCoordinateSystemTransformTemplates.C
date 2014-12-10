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

#include "fieldCoordinateSystemTransform.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Time.H"
#include "transformGeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::fieldCoordinateSystemTransform::transformField
(
    const Type& field
) const
{
    const word& fieldName = field.name() + ":Transformed";

    if (!obr_.foundObject<Type>(fieldName))
    {
        obr_.store
        (
            new Type
            (
                IOobject
                (
                    fieldName,
                    obr_.time().timeName(),
                    obr_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                field
            )
        );
    }

    Type& transField =
        const_cast<Type&>(obr_.lookupObject<Type>(fieldName));

    transField == field;

    dimensionedTensor R("R", field.dimensions(), coordSys_.R().R());

    Foam::transform(transField, R, transField);

    Info<< "    writing field " << transField.name() << nl << endl;

    transField.write();
}


template<class Type>
void Foam::fieldCoordinateSystemTransform::transform
(
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> vfType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sfType;

    if (obr_.foundObject<vfType>(fieldName))
    {
        if (debug)
        {
            Info<< type() << ": Field " << fieldName << " already in database"
                << endl;
        }

        transformField<vfType>(obr_.lookupObject<vfType>(fieldName));
    }
    else if (obr_.foundObject<sfType>(fieldName))
    {
        if (debug)
        {
            Info<< type() << ": Field " << fieldName << " already in database"
                << endl;
        }

        transformField<sfType>(obr_.lookupObject<sfType>(fieldName));
    }
    else
    {
        IOobject fieldHeader
        (
            fieldName,
            obr_.time().timeName(),
            obr_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if
        (
            fieldHeader.headerOk()
         && fieldHeader.headerClassName() == vfType::typeName
        )
        {
            if (debug)
            {
                Info<< type() << ": Field " << fieldName << " read from file"
                    << endl;
            }

            transformField<vfType>(obr_.lookupObject<vfType>(fieldName));
        }
        else if
        (
            fieldHeader.headerOk()
         && fieldHeader.headerClassName() == sfType::typeName
        )
        {
            if (debug)
            {
                Info<< type() << ": Field " << fieldName << " read from file"
                    << endl;
            }

            transformField<sfType>(obr_.lookupObject<sfType>(fieldName));
        }
    }
}


// ************************************************************************* //
