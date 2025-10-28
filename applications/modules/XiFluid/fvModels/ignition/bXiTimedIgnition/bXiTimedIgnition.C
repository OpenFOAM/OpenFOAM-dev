/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "bXiTimedIgnition.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(bXiTimedIgnition, 0);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::bXiTimedIgnition::readCoeffs(const dictionary& dict)
{
    start_.read(dict, mesh().time().userUnits());
    duration_.read(dict, mesh().time().userUnits());
    period_.readIfPresent(dict, mesh().time().userUnits());
    combustionDuration_.readIfPresent(dict, mesh().time().userUnits());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::bXiTimedIgnition::bXiTimedIgnition
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    bXiIgnition(name, modelType, mesh, dict),
    start_("start", mesh().time().userUnits(), dict),
    duration_("duration", mesh().time().userUnits(), dict),
    period_("period", mesh().time().userUnits(), dict, vGreat),
    combustionDuration_
    (
        "combustionDuration", mesh().time().userUnits(), dict, vGreat
    ),
    reset_(!ignited())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::bXiTimedIgnition::addSupFields() const
{
    if (ignited())
    {
        return wordList({"b"});
    }
    else
    {
        return wordList::null();
    }
}


Foam::scalar Foam::fv::bXiTimedIgnition::ignRelTime
(
    const scalar t
) const
{
    return std::fmod(t - start_.value(), period_.value());
}


bool Foam::fv::bXiTimedIgnition::igniting
(
    const dimensionedScalar duration
) const
{
    const scalar curTime = mesh().time().value();
    const scalar deltaT = mesh().time().deltaTValue();

    const bool igniting
    (
        ignRelTime(curTime) > -0.5*deltaT
     && ignRelTime(curTime) < max(duration.value(), 0.5*deltaT)
    );

    if (igniting)
    {
        reset_ = false;
    }

    return igniting;
}


bool Foam::fv::bXiTimedIgnition::igniting() const
{
    return igniting(duration_);
}


bool Foam::fv::bXiTimedIgnition::ignited() const
{
    const scalar curTime = mesh().time().value();
    const scalar deltaT = mesh().time().deltaTValue();

    const bool ignited
    (
        ignRelTime(curTime) > -0.5*deltaT
     && ignRelTime(curTime) < combustionDuration_.value() + 0.5*deltaT
    );

    if
    (
        !reset_
     && ignRelTime(curTime) > combustionDuration_.value() + 0.5*deltaT
    )
    {
        reset_ = true;

        mesh().lookupObjectRef<solvers::XiFluid>(solver::typeName).reset();
    }

    return ignited;
}


bool Foam::fv::bXiTimedIgnition::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }

    return false;
}


// ************************************************************************* //
