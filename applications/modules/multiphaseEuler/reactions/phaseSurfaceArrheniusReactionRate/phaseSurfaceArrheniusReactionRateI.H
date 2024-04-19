/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "phaseSurfaceArrheniusReactionRate.H"
#include "phaseSystem.H"
#include "diameterModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::phaseSurfaceArrheniusReactionRate::
phaseSurfaceArrheniusReactionRate
(
    const speciesTable& species,
    const objectRegistry& ob,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    ArrheniusReactionRate(species, dims, dict),
    phaseName_(dict.lookup("phase")),
    ob_(ob),
    tAv_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::phaseSurfaceArrheniusReactionRate::preEvaluate() const
{
    ArrheniusReactionRate::preEvaluate();

    const phaseModel& phase =
        ob_.lookupObject<phaseModel>(IOobject::groupName("alpha", phaseName_));

    tAv_ = phase.diameter().Av();
}


inline void Foam::phaseSurfaceArrheniusReactionRate::postEvaluate() const
{
    ArrheniusReactionRate::postEvaluate();

    tAv_.clear();
}


inline Foam::scalar Foam::phaseSurfaceArrheniusReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return ArrheniusReactionRate::operator()(p, T, c, li)*tAv_()[li];
}


inline Foam::scalar Foam::phaseSurfaceArrheniusReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return ArrheniusReactionRate::ddT(p, T, c, li)*tAv_()[li];
}


inline void Foam::phaseSurfaceArrheniusReactionRate::write(Ostream& os) const
{
    ArrheniusReactionRate::write(os);
    writeEntry(os, "phase", phaseName_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const phaseSurfaceArrheniusReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
