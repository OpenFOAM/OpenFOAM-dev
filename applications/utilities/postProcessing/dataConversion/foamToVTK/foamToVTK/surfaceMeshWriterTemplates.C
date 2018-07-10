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

#include "surfaceMeshWriter.H"
#include "writeFuns.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Field<Type>> Foam::surfaceMeshWriter::getFaceField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sfld
) const
{
    const polyBoundaryMesh& patches = sfld.mesh().boundaryMesh();

    tmp<Field<Type>> tfld(new Field<Type>(pp_.size()));
    Field<Type>& fld = tfld.ref();

    forAll(pp_.addressing(), i)
    {
        label facei = pp_.addressing()[i];

        label patchi = patches.whichPatch(facei);

        if (patchi == -1)
        {
            fld[i] = sfld[facei];
        }
        else
        {
            label localFacei = facei - patches[patchi].start();
            fld[i] = sfld.boundaryField()[patchi][localFacei];
        }
    }

    return tfld;
}


template<class Type>
void Foam::surfaceMeshWriter::write
(
    const UPtrList
    <
        const GeometricField<Type, fvsPatchField, surfaceMesh>
    >& sflds
)
{
    forAll(sflds, fieldi)
    {
        const GeometricField<Type, fvsPatchField, surfaceMesh>& fld =
            sflds[fieldi];

        os_ << fld.name() << ' ' << pTraits<Type>::nComponents << ' '
            << pp_.size() << " float" << std::endl;

        DynamicList<floatScalar> fField(pTraits<Type>::nComponents*pp_.size());
        writeFuns::insert(getFaceField(fld)(), fField);
        writeFuns::write(os_, binary_, fField);
    }
}


// ************************************************************************* //
