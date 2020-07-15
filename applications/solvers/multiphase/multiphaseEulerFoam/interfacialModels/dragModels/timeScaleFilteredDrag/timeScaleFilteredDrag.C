/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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

#include "timeScaleFilteredDrag.H"
#include "phasePair.H"
#include "swarmCorrection.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(timeScaleFilteredDrag, 0);
    addToRunTimeSelectionTable(dragModel, timeScaleFilteredDrag, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::timeScaleFilteredDrag::timeScaleFilteredDrag
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict.subDict("dragModel"), pair, registerObject),
    dragModel_
    (
        dragModel::New(dict.subDict("dragModel"), pair)
    ),
    minRelaxTime_("minRelaxTime", dimTime, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::timeScaleFilteredDrag::~timeScaleFilteredDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dragModels::timeScaleFilteredDrag::CdRe() const
{
    const volScalarField limit
    (
        sqr(pair_.dispersed().d())
       *pair_.dispersed().rho()
       /0.75
       /swarmCorrection_->Cs()
       /pair_.continuous().rho()
       /pair_.continuous().thermo().nu()
       /minRelaxTime_
    );

    return min(dragModel_->CdRe(), limit);
}


// ************************************************************************* //
