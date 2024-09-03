/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "multiCycleConstantbXiIgnition.H"
#include "psiuMulticomponentThermo.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(multiCycleConstantbXiIgnition, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            multiCycleConstantbXiIgnition,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::multiCycleConstantbXiIgnition::readCoeffs(const dictionary& dict)
{
    period_.read(dict, mesh().time().userUnits());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::multiCycleConstantbXiIgnition::multiCycleConstantbXiIgnition
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    constantbXiIgnition(name, modelType, mesh, dict),
    period_("period", mesh().time().userUnits(), dict),
    combustionDuration_("combustionDuration", mesh().time().userUnits(), dict),
    reset_(!ignited())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fv::multiCycleConstantbXiIgnition::ignRelTime
(
    const scalar t
) const
{
    return std::fmod(t - start_.value(), period_.value());
}


bool Foam::fv::multiCycleConstantbXiIgnition::igniting() const
{
    const scalar curTime = mesh().time().value();
    const scalar deltaT = mesh().time().deltaTValue();

    const bool igniting
    (
        ignRelTime(curTime) > -0.5*deltaT
     && ignRelTime(curTime) < max(duration_.value(), 0.5*deltaT)
    );

    if (igniting)
    {
        reset_ = false;
    }

    return igniting;
}


bool Foam::fv::multiCycleConstantbXiIgnition::ignited() const
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
     && ignRelTime(curTime) > combustionDuration_.value() - 0.5*deltaT
    )
    {
        reset_ = true;

        psiuMulticomponentThermo& thermo
        (
            mesh().lookupObjectRef<psiuMulticomponentThermo>
            (
                physicalProperties::typeName
            )
        );

        thermo.reset();
    }

    return ignited;
}


bool Foam::fv::multiCycleConstantbXiIgnition::read(const dictionary& dict)
{
    if (constantbXiIgnition::read(dict))
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
