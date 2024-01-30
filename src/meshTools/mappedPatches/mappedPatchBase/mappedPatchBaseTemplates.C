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

#include "mappedPatchBase.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mappedPatchBase::fromNeighbour(const Field<Type>& nbrFld) const
{
    if (sameUntransformedPatch())
    {
        return nbrFld;
    }

   if (nbrPatchIsMapped() && nbrMappedPatch().reMapNbr_)
    {
        treeMapPtr_.clear();
        treeNbrPatchFaceIndices_.clear();
        patchToPatchIsValid_ = false;
        nbrMappedPatch().reMapNbr_ = false;
    }

    if (usingTree_)
    {
        if (treeMapPtr_.empty())
        {
            calcMapping();
        }

        tmp<Field<Type>> tResult
        (
            new Field<Type>(nbrFld, treeNbrPatchFaceIndices_)
        );
        treeMapPtr_->distribute(tResult.ref());
        return transform_.transform().transform(tResult);
    }
    else
    {
        if
        (
            !patchToPatchIsValid_
         && !(symmetric() && nbrMappedPatch().patchToPatchIsValid_)
        )
        {
            calcMapping();
        }

        return
            transform_.transform().transform
            (
                patchToPatchIsValid_
              ? patchToPatchPtr_->tgtToSrc(nbrFld)
              : nbrMappedPatch().patchToPatchPtr_->srcToTgt(nbrFld)
            );
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mappedPatchBase::fromNeighbour(const tmp<Field<Type>>& nbrFld) const
{
    tmp<Field<Type>> tResult = fromNeighbour(nbrFld());
    nbrFld.clear();
    return tResult;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mappedPatchBase::toNeighbour(const Field<Type>& fld) const
{
    if (sameUntransformedPatch())
    {
        return fld;
    }

    if (nbrPatchIsMapped() && nbrMappedPatch().reMapNbr_)
    {
        treeMapPtr_.clear();
        treeNbrPatchFaceIndices_.clear();
        patchToPatchIsValid_ = false;
        nbrMappedPatch().reMapNbr_ = false;
    }

    if (usingTree_)
    {
        if (treeMapPtr_.empty())
        {
            calcMapping();
        }

        Field<Type> nbrFld(fld);
        treeMapPtr_->reverseDistribute(treeNbrPatchFaceIndices_.size(), nbrFld);
        tmp<Field<Type>> tResult(new Field<Type>(nbrPolyPatch().size()));
        tResult.ref().rmap(nbrFld, treeNbrPatchFaceIndices_);
        return transform_.transform().invTransform(tResult);
    }
    else
    {
        if
        (
            !patchToPatchIsValid_
         && !(symmetric() && nbrMappedPatch().patchToPatchIsValid_)
        )
        {
            calcMapping();
        }

        return
            transform_.transform().invTransform
            (
                patchToPatchIsValid_
              ? patchToPatchPtr_->srcToTgt(fld)
              : nbrMappedPatch().patchToPatchPtr_->tgtToSrc(fld)
            );
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mappedPatchBase::toNeighbour(const tmp<Field<Type>>& fld) const
{
    tmp<Field<Type>> tResult = toNeighbour(fld());
    fld.clear();
    return tResult;
}


// ************************************************************************* //
