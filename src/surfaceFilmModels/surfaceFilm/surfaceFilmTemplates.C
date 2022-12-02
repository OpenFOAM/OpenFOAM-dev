/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "surfaceFilm.H"
#include "zeroGradientFvPatchFields.H"
#include "mappedValueAndPatchInternalValueFvPatchFields.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class Type>
void Foam::surfaceFilm::toPrimary
(
    const label regionPatchi,
    Field<Type>& regionField
) const
{
    forAll(intCoupledPatchIDs_, i)
    {
        if (intCoupledPatchIDs_[i] == regionPatchi)
        {
            const mappedPatchBase& mpb =
                refCast<const mappedPatchBase>
                (
                    mesh().boundaryMesh()[regionPatchi]
                );
            regionField = mpb.reverseDistribute(regionField);
            return;
        }
    }

    FatalErrorInFunction
        << "Region patch ID " << regionPatchi << " not found in region mesh"
        << abort(FatalError);
}


template<class Type>
void Foam::surfaceFilm::toRegion
(
    const label regionPatchi,
    Field<Type>& primaryField
) const
{
    forAll(intCoupledPatchIDs_, i)
    {
        if (intCoupledPatchIDs_[i] == regionPatchi)
        {
            const mappedPatchBase& mpb =
                refCast<const mappedPatchBase>
                (
                    mesh().boundaryMesh()[regionPatchi]
                );
            primaryField = mpb.distribute(primaryField);
            return;
        }
    }

    FatalErrorInFunction
        << "Region patch ID " << regionPatchi << " not found in region mesh"
        << abort(FatalError);
}


template<class Type>
void Foam::surfaceFilm::toRegion
(
    Field<Type>& rf,
    const typename VolField<Type>::Boundary& pBf
) const
{
    forAll(intCoupledPatchIDs_, i)
    {
        const label regionPatchi = intCoupledPatchIDs_[i];
        const label primaryPatchi = primaryPatchIDs_[i];

        const polyPatch& regionPatch =
            mesh().boundaryMesh()[regionPatchi];

        const mappedPatchBase& mpb =
            refCast<const mappedPatchBase>(regionPatch);

        UIndirectList<Type>(rf, regionPatch.faceCells()) =
            mpb.distribute(pBf[primaryPatchi]);
    }
}


template<class Type>
Foam::wordList Foam::surfaceFilm::mappedFieldAndInternalPatchTypes() const
{
    wordList bTypes(mesh().boundaryMesh().size());

    bTypes = zeroGradientFvPatchField<Type>::typeName;

    forAll(intCoupledPatchIDs_, i)
    {
        bTypes[intCoupledPatchIDs_[i]] =
            mappedValueAndPatchInternalValueFvPatchField<Type>::typeName;
    }

    return bTypes;
}


// ************************************************************************* //
