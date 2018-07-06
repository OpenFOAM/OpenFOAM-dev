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

#include "rawTopoChangerFvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::rawTopoChangerFvMesh::setUnmappedValues
(
    GeometricField<Type, PatchField, GeoMesh>& fld,
    const PackedBoolList& mappedFace,
    const GeometricField<Type, PatchField, GeoMesh>& baseFld
)
{
    // Pout<< "Checking field " << fld.name() << endl;

    forAll(fld.boundaryField(), patchi)
    {
        PatchField<Type>& fvp = const_cast<PatchField<Type>&>
        (
            fld.boundaryField()[patchi]
        );

        const label start = fvp.patch().start();
        forAll(fvp, i)
        {
            if (!mappedFace[start+i])
            {
                // Pout<< "** Resetting unassigned value on patch "
                //    << fvp.patch().name()
                //    << " localface:" << i
                //    << " to:" << baseFld.boundaryField()[patchi][i] << endl;
                fvp[i] = baseFld.boundaryField()[patchi][i];
            }
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::rawTopoChangerFvMesh::zeroUnmappedValues
(
    const PackedBoolList& mappedFace
) const
{
    typedef GeometricField<Type, PatchField, GeoMesh> FieldType;

    const wordList fldNames(names(FieldType::typeName));

    forAll(fldNames, i)
    {
        // Pout<< "Checking field " << fldNames[i] << endl;

        FieldType& fld = lookupObjectRef<FieldType>(fldNames[i]);

        setUnmappedValues
        (
            fld,
            mappedFace,
            FieldType
            (
                IOobject
                (
                    "zero",
                    time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                *this,
                dimensioned<Type>("0", fld.dimensions(), Zero)
            )
        );
    }
}


// ************************************************************************* //
