/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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
    const phaseInterface& interface,
    const bool registerObject
)
:
    dispersedDragModel(dict.subDict("dragModel"), interface, registerObject),
    dragModel_
    (
        dragModel::New(dict.subDict("dragModel"), interface, false, false)
    ),
    minRelaxTime_("minRelaxTime", dimTime, dict)
{
    if (!isA<dispersedDragModel>(dragModel_()))
    {
        FatalErrorInFunction
            << "The sub-drag-model of a " << type()
            << " drag model must be for a dispersed configuration"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::timeScaleFilteredDrag::~timeScaleFilteredDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dragModels::timeScaleFilteredDrag::CdRe() const
{
    const volScalarField limit
    (
        sqr(interface_.dispersed().d())
       *interface_.dispersed().rho()
       /0.75
       /swarmCorrection_->Cs()
       /interface_.continuous().rho()
       /interface_.continuous().thermo().nu()
       /minRelaxTime_
    );

    return min(refCast<const dispersedDragModel>(dragModel_()).CdRe(), limit);
}


// ************************************************************************* //
