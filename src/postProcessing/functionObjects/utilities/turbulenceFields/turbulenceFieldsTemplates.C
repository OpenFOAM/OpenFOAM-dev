/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::turbulenceFields::processField
(
    const word& fieldName,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvalue
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    const word scopedName = modelName + ':' + fieldName;

    if (obr_.foundObject<FieldType>(scopedName))
    {
        FieldType& fld =
            const_cast<FieldType&>(obr_.lookupObject<FieldType>(scopedName));
        fld == tvalue();
    }
    else if (obr_.found(scopedName))
    {
        WarningIn
        (
            "void Foam::turbulenceFields::processField"
            "("
                "const word&, "
                "const tmp<GeometricField<Type, fvPatchField, volMesh> >&"
            ")"
        )   << "Cannot store turbulence field " << scopedName
            << " since an object with that name already exists"
            << nl << endl;
    }
    else
    {
        obr_.store
        (
            new FieldType
            (
                IOobject
                (
                    scopedName,
                    obr_.time().timeName(),
                    obr_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                tvalue
            )
        );
    }
}


// ************************************************************************* //
