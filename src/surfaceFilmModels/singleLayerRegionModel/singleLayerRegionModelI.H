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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::fvMesh&
Foam::regionModels::singleLayerRegionModel::primaryMesh() const
{
    return primaryMesh_;
}


inline const Foam::Time&
Foam::regionModels::singleLayerRegionModel::time() const
{
    return time_;
}


inline const Foam::Switch&
Foam::regionModels::singleLayerRegionModel::infoOutput() const
{
    return infoOutput_;
}


inline const Foam::word&
Foam::regionModels::singleLayerRegionModel::modelName() const
{
    return modelName_;
}


inline const Foam::fvMesh&
Foam::regionModels::singleLayerRegionModel::mesh() const
{
    return mesh_;
}


inline Foam::fvMesh& Foam::regionModels::singleLayerRegionModel::mesh()
{
    return mesh_;
}


inline const Foam::dictionary&
Foam::regionModels::singleLayerRegionModel::coeffs() const
{
    return coeffs_;
}


inline const Foam::dictionary&
Foam::regionModels::singleLayerRegionModel::solution() const
{
    return mesh().solution().dict();
}


inline const Foam::IOdictionary&
Foam::regionModels::singleLayerRegionModel::outputProperties() const
{
    if (!outputPropertiesPtr_.valid())
    {
        FatalErrorInFunction
            << "outputProperties dictionary not available"
            << abort(FatalError);
    }

    return outputPropertiesPtr_();
}


inline Foam::IOdictionary&
Foam::regionModels::singleLayerRegionModel::outputProperties()
{
    if (!outputPropertiesPtr_.valid())
    {
        FatalErrorInFunction
            << "outputProperties dictionary not available"
            << abort(FatalError);
    }

    return outputPropertiesPtr_();
}


inline bool Foam::regionModels::singleLayerRegionModel::isCoupledPatch
(
    const label regionPatchi
) const
{
    forAll(intCoupledPatchIDs_, i)
    {
        if (intCoupledPatchIDs_[i] == regionPatchi)
        {
            return true;
        }
    }

    return false;
}


inline bool Foam::regionModels::singleLayerRegionModel::isRegionPatch
(
    const label primaryPatchi
) const
{
    forAll(primaryPatchIDs_, i)
    {
        if (primaryPatchIDs_[i] == primaryPatchi)
        {
            return true;
        }
    }

    return false;
}


inline const Foam::labelList&
Foam::regionModels::singleLayerRegionModel::primaryPatchIDs() const
{
    return primaryPatchIDs_;
}


inline const Foam::labelList&
Foam::regionModels::singleLayerRegionModel::intCoupledPatchIDs() const
{
    return intCoupledPatchIDs_;
}


inline Foam::label Foam::regionModels::singleLayerRegionModel::regionPatchID
(
    const label primaryPatchID
) const
{
    forAll(primaryPatchIDs_, i)
    {
        if (primaryPatchIDs_[i] == primaryPatchID)
        {
            return intCoupledPatchIDs_[i];
        }
    }

    return -1;
}


// ************************************************************************* //
