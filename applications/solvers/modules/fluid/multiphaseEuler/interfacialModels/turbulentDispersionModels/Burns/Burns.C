/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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

#include "Burns.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "dispersedDragModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentDispersionModels
{
    defineTypeNameAndDebug(Burns, 0);
    addToRunTimeSelectionTable
    (
        turbulentDispersionModel,
        Burns,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::Burns::Burns
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    dispersedTurbulentDispersionModel(dict, interface),
    sigma_("sigma", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::Burns::~Burns()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::turbulentDispersionModels::Burns::D() const
{
    const dragModels::dispersedDragModel& drag =
        interface_.mesh().lookupObject<dragModels::dispersedDragModel>
        (
            IOobject::groupName
            (
                dragModel::typeName,
                interface_.name()
            )
        );

    return
        drag.Ki()
       *continuousTurbulence().nut()
       /sigma_
       *interface_.dispersed()
       *sqr(interface_.dispersed() + interface_.continuous())
       /(
            max
            (
                interface_.dispersed(),
                interface_.dispersed().residualAlpha()
            )
           *max
            (
                interface_.continuous(),
                interface_.continuous().residualAlpha()
            )
        );
}


// ************************************************************************* //
