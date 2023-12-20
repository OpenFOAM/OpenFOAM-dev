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

#include "patchWriter.H"
#include "vtkWriteFieldOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::patchWriter::write
(
    const UPtrList<const VolField<Type>>& flds
)
{
    forAll(flds, fieldi)
    {
        const VolField<Type>& fld = flds[fieldi];

        os_ << fld.name() << ' ' << pTraits<Type>::nComponents << ' '
            << nFaces_ << " float" << std::endl;

        DynamicList<floatScalar> fField(pTraits<Type>::nComponents*nFaces_);

        forAll(patchIndices_, j)
        {
            label patchi = patchIndices_[j];

            const fvPatchField<Type>& pfld = fld.boundaryField()[patchi];

            if (nearCellValue_)
            {
                vtkWriteOps::insert(pfld.patchInternalField()(), fField);
            }
            else
            {
                vtkWriteOps::insert(pfld, fField);
            }
        }
        vtkWriteOps::write(os_, binary_, fField);
    }
}


template<class Type>
void Foam::patchWriter::write
(
    const UPtrList<const PointField<Type>>& flds
)
{
    forAll(flds, fieldi)
    {
        const PointField<Type>& fld =
            flds[fieldi];

        os_ << fld.name() << ' ' << pTraits<Type>::nComponents << ' '
            << nPoints_ << " float" << std::endl;

        DynamicList<floatScalar> fField(pTraits<Type>::nComponents*nPoints_);

        forAll(patchIndices_, j)
        {
            label patchi = patchIndices_[j];

            const pointPatchField<Type>& pfld = fld.boundaryField()[patchi];

            vtkWriteOps::insert(pfld.patchInternalField()(), fField);
        }
        vtkWriteOps::write(os_, binary_, fField);
    }
}


template<class Type>
void Foam::patchWriter::write
(
    const PrimitivePatchInterpolation<primitivePatch>& pInter,
    const UPtrList<const VolField<Type>>& flds
)
{
    forAll(flds, fieldi)
    {
        const VolField<Type>& fld = flds[fieldi];

        os_ << fld.name() << ' ' << pTraits<Type>::nComponents << ' '
            << nPoints_ << " float" << std::endl;

        DynamicList<floatScalar> fField(pTraits<Type>::nComponents*nPoints_);

        forAll(patchIndices_, j)
        {
            label patchi = patchIndices_[j];

            const fvPatchField<Type>& pfld = fld.boundaryField()[patchi];

            if (nearCellValue_)
            {
                vtkWriteOps::insert
                (
                    pInter.faceToPointInterpolate
                    (
                        pfld.patchInternalField()()
                    )(),
                    fField
                );
            }
            else
            {
                vtkWriteOps::insert
                (
                    pInter.faceToPointInterpolate(pfld)(),
                    fField
                );
            }
        }
        vtkWriteOps::write(os_, binary_, fField);
    }
}


// ************************************************************************* //
