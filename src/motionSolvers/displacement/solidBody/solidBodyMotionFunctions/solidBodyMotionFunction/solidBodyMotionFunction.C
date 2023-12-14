/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidBodyMotionFunction, 0);
    defineRunTimeSelectionTable(solidBodyMotionFunction, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunction::solidBodyMotionFunction
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    SBMFCoeffs_
    (
        SBMFCoeffs.optionalSubDict
        (
            word(SBMFCoeffs.lookup("solidBodyMotionFunction")) + "Coeffs"
        )
    ),
    time_(runTime)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunction::~solidBodyMotionFunction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::solidBodyMotionFunction::read(const dictionary& SBMFCoeffs)
{
    SBMFCoeffs_ = SBMFCoeffs.optionalSubDict(type() + "Coeffs");

    return true;
}


void Foam::solidBodyMotionFunction::writeData(Ostream& os) const
{
    os << SBMFCoeffs_;
}


// ************************************************************************* //
