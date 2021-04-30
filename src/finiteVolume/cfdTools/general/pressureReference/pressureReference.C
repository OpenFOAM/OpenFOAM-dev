/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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

#include "pressureReference.H"
#include "findRefCell.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureReference::pressureReference
(
    const volScalarField& p,
    const volScalarField& pRef,
    const dictionary& dict,
    const bool pRefRequired
)
:
    refCell_(-1),
    refValue_(0)
{
    // Set the reference cell and value for closed domain simulations
    if (pRefRequired)
    {
        setRefCell(p, pRef, dict, refCell_, refValue_);
    }

    // Print update message fvSolution -> fvConstraints
    if
    (
        dict.found("pMin")
     || dict.found("pMax")
     || dict.found("pMinFactor")
     || dict.found("pMaxFactor")
     || dict.found("rhoMin")
     || dict.found("rhoMax")
    )
    {
        FatalIOErrorInFunction(dict)
            << "Pressure limits should now be specified in fvConstraints:\n\n"
               "limitp\n"
               "{\n"
               "    type       limitPressure;\n"
               "\n";

        if (dict.found("pMin"))
        {
            FatalIOError
                << "    min        " << dict.lookup<scalar>("pMin")
                << ";\n";
        }

        if (dict.found("pMax"))
        {
            FatalIOError
                << "    max        " << dict.lookup<scalar>("pMax")
                << ";\n";
        }

        if (dict.found("pMinFactor"))
        {
            FatalIOError
                << "    minFactor  " << dict.lookup<scalar>("pMinFactor")
                << ";\n";
        }

        if (dict.found("pMaxFactor"))
        {
            FatalIOError
                << "    maxFactor  " << dict.lookup<scalar>("pMaxFactor")
                << ";\n";
        }

        FatalIOError << "}\n";

        if
        (
            dict.found("rhoMin")
         || dict.found("rhoMax")
        )
        {
            FatalIOError
                << "\nrhoMin and rhoMax are no longer supported.\n";
        }

        FatalIOError << exit(FatalIOError);
    }
}


Foam::pressureReference::pressureReference
(
    const volScalarField& p,
    const dictionary& dict,
    const bool pRefRequired
)
:
    pressureReference(p, p, dict, pRefRequired)
{}


// ************************************************************************* //
