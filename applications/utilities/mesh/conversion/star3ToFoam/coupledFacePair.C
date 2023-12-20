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

Description
    Data associated with a pair of coupled faces.
    1 represents integral match; all other number are arbitrary matches

\*---------------------------------------------------------------------------*/

#include "coupledFacePair.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledFacePair::coupledFacePair
(
    const label coupleNo,
    const label mC, const label mF,
    const label sC, const label sF,
    const label integral
)
:
    coupleIndex_(coupleNo),
    masterCellIndex_(mC),
    masterFaceIndex_(mF),
    slaveCellIndex_(sC),
    slaveFaceIndex_(sF),
    integralMatch_(integral == 1)
{}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const coupledFacePair& c)
{
    os  << "Master cell: " << c.masterCellIndex_
        << " face: " << c.masterFaceIndex_ << endl
        << "Slave cell: " << c.slaveCellIndex_
        << " face: " << c.slaveFaceIndex_ << endl
        << "Integral: " << c.integralMatch_ << endl;

    return os;
}


// ************************************************************************* //
