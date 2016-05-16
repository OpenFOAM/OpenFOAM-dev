/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"

template<class FieldType>
FieldType& Foam::functionObjects::mag::magField
(
    const word& magName,
    const dimensionSet& dims
)
{
    if (!mesh_.foundObject<FieldType>(magName))
    {
        FieldType* magFieldPtr
        (
            new FieldType
            (
                IOobject
                (
                    magName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dims, 0.0)
            )
        );

        mesh_.objectRegistry::store(magFieldPtr);
    }

    const FieldType& f = mesh_.lookupObject<FieldType>(magName);

    return const_cast<FieldType&>(f);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::mag::calc
(
    const word& fieldName,
    const word& resultName,
    bool& processed
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> vfType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sfType;

    if (mesh_.foundObject<vfType>(fieldName))
    {
        const vfType& vf = mesh_.lookupObject<vfType>(fieldName);

        volScalarField& field =
            magField<volScalarField>(resultName_, vf.dimensions());

        field = Foam::mag(vf);

        processed = true;
    }
    else if (mesh_.foundObject<sfType>(fieldName))
    {
        const sfType& sf = mesh_.lookupObject<sfType>(fieldName);

        surfaceScalarField& field =
            magField<surfaceScalarField>(resultName_, sf.dimensions());

        field = Foam::mag(sf);

        processed = true;
    }
}


// ************************************************************************* //
