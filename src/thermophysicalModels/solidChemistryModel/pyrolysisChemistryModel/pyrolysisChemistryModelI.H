/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "volFields.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class SolidThermo, class GasThermo>
inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::pyrolysisChemistryModel<CompType, SolidThermo, GasThermo>::RRg()
{
    return RRg_;
}


template<class CompType, class SolidThermo, class GasThermo>
inline const Foam::PtrList<GasThermo>&
Foam::pyrolysisChemistryModel<CompType, SolidThermo, GasThermo>::
gasThermo() const
{
    return gasThermo_;
}


template<class CompType, class SolidThermo, class GasThermo>
inline const Foam::speciesTable&
Foam::pyrolysisChemistryModel<CompType, SolidThermo, GasThermo>::
gasTable() const
{
    return pyrolisisGases_;
}


template<class CompType, class SolidThermo, class GasThermo>
inline Foam::label
Foam::pyrolysisChemistryModel<CompType, SolidThermo, GasThermo>::
nSpecie() const
{
    return nSpecie_;
}


template<class CompType, class SolidThermo, class GasThermo>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::pyrolysisChemistryModel<CompType, SolidThermo, GasThermo>::RRg
(
    const label i
) const
{
    return RRg_[i];
}


template<class CompType, class SolidThermo, class GasThermo>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::pyrolysisChemistryModel<CompType, SolidThermo, GasThermo>::
RRg() const
{
    tmp<volScalarField::Internal> tRRg
    (
        volScalarField::Internal::New
        (
            "RRg",
            this->mesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, 0)
        )
    );

    if (this->chemistry_)
    {
        volScalarField::Internal& RRg = tRRg.ref();
        for (label i=0; i < nGases_; i++)
        {
            RRg += RRg_[i];
        }
    }

    return tRRg;
}


// ************************************************************************* //
