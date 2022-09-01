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

#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mappedPatchBase::distribute(const Field<Type>& fld) const
{
    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            if (AMIPtr_.empty())
            {
                calcAMI();
            }
            return AMIPtr_->interpolateToSource(fld);
        }
        case PATCHTOPATCH:
        {
            if
            (
                !patchToPatchIsValid_
             && !(
                    sampleIsMappedPatch()
                 && sampleMappedPatch().patchToPatchIsValid_
                )
            )
            {
                calcPatchToPatch();
            }

            return
                patchToPatchIsValid_
              ? patchToPatchPtr_->tgtToSrc(fld)
              : sampleMappedPatch().patchToPatchPtr_->srcToTgt(fld);
        }
        default:
        {
            if (mapPtr_.empty())
            {
                calcMapping();
            }
            tmp<Field<Type>> tResult(new Field<Type>(fld, mapIndices_));
            mapPtr_->distribute(tResult.ref());
            return tResult;
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mappedPatchBase::distribute(const tmp<Field<Type>>& fld) const
{
    tmp<Field<Type>> tResult = distribute(fld());
    fld.clear();
    return tResult;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mappedPatchBase::reverseDistribute(const Field<Type>& fld) const
{
    switch (mode_)
    {
        case NEARESTPATCHFACE:
        {
            if (mapPtr_.empty())
            {
                calcMapping();
            }
            Field<Type> sampleFld(fld);
            mapPtr_->reverseDistribute(mapIndices_.size(), sampleFld);
            tmp<Field<Type>> tResult(new Field<Type>(samplePolyPatch().size()));
            UIndirectList<Type>(tResult.ref(), mapIndices_) = fld;
            return tResult;
        }
        case NEARESTPATCHFACEAMI:
        {
            if (AMIPtr_.empty())
            {
                calcAMI();
            }
            return AMIPtr_->interpolateToTarget(fld);
        }
        case PATCHTOPATCH:
        {
            if
            (
                !patchToPatchIsValid_
             && !(
                    sampleIsMappedPatch()
                 && sampleMappedPatch().patchToPatchIsValid_
                )
            )
            {
                calcPatchToPatch();
            }

            return
                patchToPatchIsValid_
              ? patchToPatchPtr_->srcToTgt(fld)
              : sampleMappedPatch().patchToPatchPtr_->tgtToSrc(fld);

        }
        default:
        {
            FatalErrorInFunction
                << "Reverse distribute can only be used in "
                << sampleModeNames_[NEARESTPATCHFACE] << ", "
                << sampleModeNames_[NEARESTPATCHFACEAMI] << " or "
                << sampleModeNames_[PATCHTOPATCH] << " mode"
                << exit(FatalError);

            return tmp<Field<Type>>(nullptr);
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mappedPatchBase::reverseDistribute(const tmp<Field<Type>>& fld) const
{
    tmp<Field<Type>> tResult = reverseDistribute(fld());
    fld.clear();
    return tResult;
}


// ************************************************************************* //
