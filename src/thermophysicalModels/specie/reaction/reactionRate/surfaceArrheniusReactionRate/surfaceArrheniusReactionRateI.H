/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2024 OpenFOAM Foundation
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

#include "surfaceArrheniusReactionRate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::surfaceArrheniusReactionRate::surfaceArrheniusReactionRate
(
    const speciesTable& species,
    const objectRegistry& ob,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    ArrheniusReactionRate(species, dims, dict),
    AvName_(dict.lookupBackwardsCompatible({"Av", "a"})),
    ob_(ob),
    tAv_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::surfaceArrheniusReactionRate::preEvaluate() const
{
    ArrheniusReactionRate::preEvaluate();

    tAv_ = ob_.lookupObject<volScalarField::Internal>(AvName_);
}


inline void Foam::surfaceArrheniusReactionRate::postEvaluate() const
{
    ArrheniusReactionRate::postEvaluate();

    tAv_ = tmp<volScalarField::Internal>(nullptr);
}


inline Foam::scalar Foam::surfaceArrheniusReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return ArrheniusReactionRate::operator()(p, T, c, li)*tAv_()[li];
}


inline Foam::scalar Foam::surfaceArrheniusReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return ArrheniusReactionRate::ddT(p, T, c, li)*tAv_()[li];
}


inline void Foam::surfaceArrheniusReactionRate::write(Ostream& os) const
{
    ArrheniusReactionRate::write(os);
    writeEntry(os, "a", AvName_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const surfaceArrheniusReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
