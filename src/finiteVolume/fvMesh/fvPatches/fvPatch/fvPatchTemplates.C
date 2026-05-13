/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "fvPatch.H"
#include "SlicedDimensionedField.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvPatch::patchInternalField
(
    const UList<Type>& f
) const
{
    tmp<Field<Type>> tpif(new Field<Type>(size()));
    Field<Type>& pif = tpif.ref();

    const labelUList& faceCells = this->faceCells();

    forAll(pif, facei)
    {
        pif[facei] = f[faceCells[facei]];
    }

    return tpif;
}


template<class Type>
void Foam::fvPatch::patchInternalField
(
    const UList<Type>& f,
    Field<Type>& pif
) const
{
    pif.setSize(size());

    const labelUList& faceCells = this->faceCells();

    forAll(pif, facei)
    {
        pif[facei] = f[faceCells[facei]];
    }
}


template<class GeometricField, class Type>
const typename GeometricField::Patch& Foam::fvPatch::patchField
(
    const GeometricField& gf
) const
{
    return gf.boundaryField()[index()];
}


template<class GeometricField, class Type>
typename GeometricField::Patch& Foam::fvPatch::patchField
(
    GeometricField& gf
) const
{
    return gf.boundaryFieldRef()[index()];
}


template<class GeometricField, class Type>
const typename GeometricField::Patch& Foam::fvPatch::lookupPatchField
(
    const word& name
) const
{
    return patchField<GeometricField, Type>
    (
        db().template lookupObject<GeometricField>(name)
    );
}


template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::fvPatch>>
Foam::fvPatch::lookupField
(
    const word& name
) const
{
    typedef GeometricField<Type, fvMesh> volFieldType;
    typedef GeometricField<Type, surfaceMesh> surfaceFieldType;

    const bool haveVf = db().template foundObject<volFieldType>(name);
    const bool haveSf = db().template foundObject<surfaceFieldType>(name);

    if (!haveVf && !haveSf)
    {
        FatalErrorInFunction
            << nl << "    request for " << volFieldType::typeName
            << " or " << surfaceFieldType::typeName << " " << name
            << " from objectRegistry " << db().name()
            << " failed\n    available objects of type "
            << volFieldType::typeName << " are" << nl
            << db().toc<volFieldType>()
            << "\n    available objects of type "
            << surfaceFieldType::typeName << " are" << nl
            << db().toc<surfaceFieldType>()
            << abort(FatalError);
    }

    const IOobject io
    (
        name + '_' + this->name(),
        time().name(),
        db(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );

    if (haveVf)
    {
        const volFieldType& vf =
            db().template lookupObject<volFieldType>(name);

        return tmp<DimensionedField<Type, fvPatch>>
        (
            new SlicedDimensionedField<Type, fvPatch>
            (
                io,
                *this,
                vf.dimensions(),
                vf.boundaryField()[index()]
            ),
            true
        );
    }
    else // if (haveSf)
    {
        const surfaceFieldType& sf =
            db().template lookupObject<surfaceFieldType>(name);

        return tmp<DimensionedField<Type, fvPatch>>
        (
            new SlicedDimensionedField<Type, fvPatch>
            (
                io,
                *this,
                sf.dimensions(),
                sf.boundaryField()[index()]
            ),
            true
        );
    }
}


// ************************************************************************* //
