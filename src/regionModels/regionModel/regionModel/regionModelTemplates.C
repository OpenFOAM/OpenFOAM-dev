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

#include "regionModel.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::regionModels::regionModel::toPrimary
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
                    regionMesh().boundaryMesh()[regionPatchi]
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
void Foam::regionModels::regionModel::toRegion
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
                    regionMesh().boundaryMesh()[regionPatchi]
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
void Foam::regionModels::regionModel::toRegion
(
    Field<Type>& rf,
    const typename GeometricField<Type, fvPatchField, volMesh>::Boundary& pBf
) const
{
    forAll(intCoupledPatchIDs_, i)
    {
        const label regionPatchi = intCoupledPatchIDs_[i];
        const label primaryPatchi = primaryPatchIDs_[i];

        const polyPatch& regionPatch =
            regionMesh().boundaryMesh()[regionPatchi];

        const mappedPatchBase& mpb =
            refCast<const mappedPatchBase>(regionPatch);

        UIndirectList<Type>(rf, regionPatch.faceCells()) =
            mpb.distribute(pBf[primaryPatchi]);
    }
}


// ************************************************************************* //
