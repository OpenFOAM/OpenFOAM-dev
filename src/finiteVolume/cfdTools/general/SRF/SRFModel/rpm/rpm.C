/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "rpm.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace SRF
    {
        defineTypeNameAndDebug(rpm, 0);

        addToRunTimeSelectionTable
        (
            SRFModel,
            rpm,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SRF::rpm::rpm
(
    const volVectorField& U
)
:
    SRFModel(typeName, U),
    rpm_(readScalar(SRFModelCoeffs_.lookup("rpm")))
{
    // Initialise the angular velocity
    omega_.value() = axis_*rpm_*constant::mathematical::twoPi/60.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SRF::rpm::~rpm()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::SRF::rpm::read()
{
    if (SRFModel::read())
    {
        // Re-read rpm
        SRFModelCoeffs_.lookup("rpm") >> rpm_;

        // Update angular velocity
        omega_.value() = axis_*rpm_*(constant::mathematical::twoPi/60.0);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
