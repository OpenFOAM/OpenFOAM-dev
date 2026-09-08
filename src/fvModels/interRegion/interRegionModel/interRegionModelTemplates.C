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

#include "interRegionModel.H"
#include "fvModels.H"

// * * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * //

template<class InterRegionModelType>
const InterRegionModelType& Foam::fv::interRegionModel::nbrModel() const
{
    const Foam::fvModels& nbrModels = fvModels::New(nbrMesh());

    wordList nbrModelNames;
    label nbrModeli = -1;

    forAll(nbrModels, i)
    {
        if (isA<InterRegionModelType>(nbrModels[i]))
        {
            const InterRegionModelType& model =
                refCast<const InterRegionModelType>(nbrModels[i]);

            nbrModelNames.append(model.name());

            if
            (
                model.nbrRegionName() == mesh().name()
             && (nbrModelName_.empty() || model.name() == nbrModelName_)
            )
            {
                nbrModeli = nbrModeli == -1 ? i : -2;
            }
        }
    }

    if (nbrModeli < 0 && !nbrModelNames.empty() && !nbrModelName_.empty())
    {
        FatalErrorInFunction
            << "Neighbour for model '" << name() << "' in region '"
            << mesh().name() << "' not found. No model of type '"
            << InterRegionModelType::typeName << "' named '" << nbrModelName_
            << "' was found in region '" << nbrMesh().name() << "'." << nl << nl
            << "Available models of type '" << InterRegionModelType::typeName
            << "' in region '" << nbrMesh().name() << "' are:" << nl
            << nbrModelNames << exit(FatalError);
    }

    if (nbrModeli == -1)
    {
        FatalErrorInFunction
            << "Neighbour for model '" << name() << "' in region '"
            << mesh().name() << "' not found. No models of type '"
            << InterRegionModelType::typeName << "' were found in neighbour "
            << "region '" << nbrMesh().name() << "'." << exit(FatalError);
    }

    if (nbrModeli == -2)
    {
        nbrModelNameDict_.lookup<word>("neighbourModel");
    }

    const InterRegionModelType& nbrModel =
        refCast<const InterRegionModelType>(nbrModels[nbrModeli]);

    if
    (
        nbrModel.nbrRegionName() != mesh().name()
     || (!nbrModel.nbrModelName().empty() && nbrModel.nbrModelName() != name())
    )
    {
        FatalErrorInFunction
            << "Model '" << name() << "' in region '" << mesh().name()
            << "' neighbours model '" << nbrModel.name() << "' in region '"
            << nbrMesh().name() << "' but the reverse is not true"
            << exit(FatalError);
    }

    // Check and synchronise the ownership
    if (owner_ == nbrModel.owner_)
    {
        FatalErrorInFunction
            << (owner_ == +1 ? "Both" : "Neither") << " of corresponding "
            << "models '" << name() << "' in region '" << mesh().name()
            << "' and '" << nbrModel.name() << "' in region '"
            << nbrMesh().name() << "' "
            << (owner_ == -1 ? "specify" : "are specified to be")
            << " the owner" << exit(FatalError);
    }
    if (owner_ == -1)
    {
        owner_ = 1 - nbrModel.owner_;
    }
    if (nbrModel.owner_ == -1)
    {
        nbrModel.owner_ = 1 - owner_;
    }

    return refCast<const InterRegionModelType>(nbrModels[nbrModeli]);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fv::interRegionModel::interpolate
(
    const Field<Type>& field
) const
{
    if (!neighbour())
    {
        return interpolation().tgtToSrc(field);
    }
    else
    {
        return nbrModel().interpolation().srcToTgt(field);
    }
}


template<class Type>
void Foam::fv::interRegionModel::interpolate
(
    const Field<Type>& field,
    Field<Type>& result
) const
{
    if (!neighbour())
    {
        result = interpolation().tgtToSrc(field, result);
    }
    else
    {
        result = nbrModel().interpolation().srcToTgt(field, result);
    }
}


// ************************************************************************* //
