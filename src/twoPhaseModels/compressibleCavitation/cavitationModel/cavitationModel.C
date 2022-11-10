/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "cavitationModel.H"
#include "constantPressure.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
    defineTypeNameAndDebug(cavitationModel, 0);
    defineRunTimeSelectionTable(cavitationModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressible::cavitationModel::cavitationModel
(
    const dictionary& dict,
    const compressibleTwoPhases& phases
)
:
    phases_(phases),
    liquidIndex_(phases.index(dict.lookup<word>("liquid"))),
    saturationModel_(saturationPressureModel::New("pSat", dict))
{}


bool Foam::compressible::cavitationModel::read(const dictionary& dict)
{
    saturationModel_.reset(saturationPressureModel::New("pSat", dict).ptr());

    return true;
}


// ************************************************************************* //
