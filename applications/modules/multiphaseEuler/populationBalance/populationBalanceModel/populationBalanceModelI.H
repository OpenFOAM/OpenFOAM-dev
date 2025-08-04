/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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

#include "populationBalanceModel.H"
#include "populationBalance.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

inline Foam::label Foam::populationBalanceModel::sourceUpdateInterval() const
{
    return solverDict().lookupOrDefault<label>("sourceUpdateInterval", 1);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::phaseSystem& Foam::populationBalanceModel::fluid() const
{
    return fluid_;
}


inline const Foam::fvMesh& Foam::populationBalanceModel::mesh() const
{
    return mesh_;
}


inline const Foam::phaseModel&
Foam::populationBalanceModel::continuousPhase() const
{
    return continuousPhase_;
}


inline Foam::label Foam::populationBalanceModel::nGroups() const
{
    return phases_.size();
}


inline const Foam::UPtrList<const Foam::phaseModel>&
Foam::populationBalanceModel::phases() const
{
    return phases_;
}


inline const Foam::UPtrList<const Foam::phaseModel>&
Foam::populationBalanceModel::uniquePhases() const
{
    return uniquePhases_;
}


inline const Foam::UPtrList<const Foam::diameterModels::populationBalance>&
Foam::populationBalanceModel::diameters() const
{
    return diameters_;
}


inline const Foam::UPtrList<const Foam::diameterModels::populationBalance>&
Foam::populationBalanceModel::uniqueDiameters() const
{
    return uniqueDiameters_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::populationBalanceModel::fs() const
{
    return fs_;
}


inline Foam::PtrList<Foam::volScalarField>& Foam::populationBalanceModel::fs()
{
    return fs_;
}


inline const Foam::volScalarField& Foam::populationBalanceModel::f
(
    const label i
) const
{
    return fs_[i];
}


inline Foam::volScalarField& Foam::populationBalanceModel::f(const label i)
{
    return fs_[i];
}


inline const Foam::PtrList<Foam::dimensionedScalar>&
Foam::populationBalanceModel::dSphs() const
{
    return dSphs_;
}


inline const Foam::dimensionedScalar& Foam::populationBalanceModel::dSph
(
    const label i
) const
{
    return dSphs_[i];
}


inline const Foam::PtrList<Foam::dimensionedScalar>&
Foam::populationBalanceModel::vs() const
{
    return vs_;
}


inline const Foam::dimensionedScalar& Foam::populationBalanceModel::v
(
    const label i
) const
{
    return vs_[i];
}


inline const Foam::populationBalanceModel::dmdtfTable&
Foam::populationBalanceModel::dmdtfs() const
{
    return dmdtfs_;
}


inline const Foam::populationBalanceModel::dmdtfTable&
Foam::populationBalanceModel::expansionDmdtfs() const
{
    return expansionDmdtfs_;
}


inline const Foam::populationBalanceModel::dmdtfTable&
Foam::populationBalanceModel::modelSourceDmdtfs() const
{
    return modelSourceDmdtfs_;
}


inline const Foam::dictionary& Foam::populationBalanceModel::solverDict() const
{
    return mesh_.solution().solverDict(name_);
}


inline bool Foam::populationBalanceModel::solveOnFinalIterOnly() const
{
    return solverDict().lookupOrDefault<bool>("solveOnFinalIterOnly", false);
}


inline const Foam::List<Foam::Pair<Foam::label>>&
Foam::populationBalanceModel::coalescencePairs() const
{
    return coalescencePairs_;
}


inline const Foam::List<Foam::Pair<Foam::label>>&
Foam::populationBalanceModel::binaryBreakupPairs() const
{
    return binaryBreakupPairs_;
}


inline const Foam::volScalarField& Foam::populationBalanceModel::alphas() const
{
    return
        uniquePhases_.size() > 1
      ? alphas_()
      : uniquePhases_.first();
}


inline Foam::tmp<Foam::volVectorField> Foam::populationBalanceModel::U() const
{
    return
        uniquePhases_.size() > 1
      ? tmp<volVectorField>(U_())
      : uniquePhases_.first().U();
}


// ************************************************************************* //
