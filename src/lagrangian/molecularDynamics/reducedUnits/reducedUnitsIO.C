/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "reducedUnits.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const reducedUnits& rU)
{
    os  << nl << "Defined: " << nl
        << tab << "refLength = " << rU.refLength() << " m" << nl
        << tab << "refTime = " << rU.refTime() << " s" << nl
        << tab << "refMass = " << rU.refMass() << " kg" << nl
        << tab << "Boltzmann constant, kb = " << reducedUnits::kb << " J/K"
        << nl << "Calculated: " << nl
        << tab << "refEnergy = " << rU.refEnergy() << " J" << nl
        << tab << "refTemp = " << rU.refTemp() << " K" << nl
        << tab << "refForce = " << rU.refForce() << " N" << nl
        << tab << "refVelocity = " << rU.refVelocity() << " m/s" << nl
        << tab << "refVolume = " << rU.refVolume() << " m^3" << nl
        << tab << "refPressure = " << rU.refPressure() << " N/m^2" << nl
        << tab << "refMassDensity = " << rU.refMassDensity() << " kg/m^3" << nl
        << tab << "refNumberDensity = " << rU.refNumberDensity() << " m^-3"
        << endl;

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<(Foam::Ostream&, "
        "const Foam::reducedUnits&)"
    );

    return os;
}


// ************************************************************************* //
