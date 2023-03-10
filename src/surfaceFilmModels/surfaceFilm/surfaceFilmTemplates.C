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

#include "surfaceFilm.H"
#include "zeroGradientFvPatchFields.H"
#include "mappedValueAndPatchInternalValueFvPatchFields.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class Type>
void Foam::surfaceFilm::toPrimary
(
    const label filmPatchi,
    Field<Type>& filmField
) const
{
    forAll(intCoupledPatchIDs_, i)
    {
        if (intCoupledPatchIDs_[i] == filmPatchi)
        {
            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>
                (
                    mesh().boundaryMesh()[filmPatchi]
                );
            filmField = mpp.toNeigbour(filmField);
            return;
        }
    }

    FatalErrorInFunction
        << "Film patch ID " << filmPatchi << " not found in film mesh"
        << abort(FatalError);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::surfaceFilm::toFilm
(
    const label filmPatchi,
    const Field<Type>& primaryPatchField
) const
{
    forAll(intCoupledPatchIDs_, i)
    {
        if (intCoupledPatchIDs_[i] == filmPatchi)
        {
            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>
                (
                    mesh().boundaryMesh()[filmPatchi]
                );

            return mpp.fromNeigbour(primaryPatchField);
        }
    }

    FatalErrorInFunction
        << "Film patch ID " << filmPatchi << " not found in film mesh"
        << abort(FatalError);

    return tmp<Field<Type>>(nullptr);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::surfaceFilm::toFilm
(
    const label filmPatchi,
    const tmp<Field<Type>>& tprimaryPatchField
) const
{
    forAll(intCoupledPatchIDs_, i)
    {
        if (intCoupledPatchIDs_[i] == filmPatchi)
        {
            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>
                (
                    mesh().boundaryMesh()[filmPatchi]
                );

            return mpp.fromNeigbour(tprimaryPatchField);
        }
    }

    FatalErrorInFunction
        << "Film patch ID " << filmPatchi << " not found in film mesh"
        << abort(FatalError);

    return tmp<Field<Type>>(nullptr);
}


template<class Type>
void Foam::surfaceFilm::toFilm
(
    Field<Type>& rf,
    const typename VolField<Type>::Boundary& pBf
) const
{
    forAll(intCoupledPatchIDs_, i)
    {
        const label filmPatchi = intCoupledPatchIDs_[i];
        const label primaryPatchi = primaryPatchIDs_[i];

        const polyPatch& filmPatch =
            mesh().boundaryMesh()[filmPatchi];

        const mappedPatchBase& mpp =
            refCast<const mappedPatchBase>(filmPatch);

        UIndirectList<Type>(rf, filmPatch.faceCells()) =
            mpp.fromNeigbour(pBf[primaryPatchi]);
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
