/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Copyright (C) 2022 VTT Technical Research Centre of Finland Ltd
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "residualDiameter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(residualDiameter, 0);

    addToRunTimeSelectionTable
    (
        diameterModel,
        residualDiameter,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::residualDiameter::residualDiameter
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    spherical(diameterProperties, phase),
    d_("d", dimLength, diameterProperties),
    dResidual_("dResidual", dimLength, diameterProperties)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::residualDiameter::~residualDiameter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::diameterModels::residualDiameter::
d() const
{
    return
        pos0(phase()-phase().residualAlpha())*d_
      + neg(phase()-phase().residualAlpha())*dResidual_;
}


bool Foam::diameterModels::residualDiameter::read
(
    const dictionary& phaseProperties
)
{
    spherical::read(phaseProperties);

    diameterProperties().lookup("d") >> d_;
    diameterProperties().lookup("dResidual") >> dResidual_;

    return true;
}


// ************************************************************************* //
