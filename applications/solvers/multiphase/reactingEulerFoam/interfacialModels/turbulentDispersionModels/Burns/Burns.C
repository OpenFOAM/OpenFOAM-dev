/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2020 OpenFOAM Foundation
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
#include "phasePair.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

#include "dragModel.H"

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
    const phasePair& pair
)
:
    turbulentDispersionModel(dict, pair),
    sigma_("sigma", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::Burns::~Burns()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::turbulentDispersionModels::Burns::D() const
{
    const fvMesh& mesh(pair_.phase1().mesh());
    const dragModel& drag =
        mesh.lookupObject<dragModel>
        (
            IOobject::groupName(dragModel::typeName, pair_.name())
        );

    return
        drag.Ki()
       *continuousTurbulence().nut()
       /sigma_
       *pair_.dispersed()
       *sqr(pair_.dispersed() + pair_.continuous())
       /(
            max(pair_.dispersed(), pair_.dispersed().residualAlpha())
           *max(pair_.continuous(), pair_.continuous().residualAlpha())
        );
}


// ************************************************************************* //
