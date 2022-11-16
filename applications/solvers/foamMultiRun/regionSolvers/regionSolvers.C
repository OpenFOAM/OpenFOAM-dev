/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "regionSolvers.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolvers::regionSolvers(const Time& runTime)
{
    if (runTime.controlDict().found("regionSolvers"))
    {
        const dictionary& regionSolversDict =
            runTime.controlDict().subDict("regionSolvers");

        forAllConstIter(dictionary, regionSolversDict, iter)
        {
            append(Pair<word>(iter().keyword(), iter().stream()));
        }
    }
    else
    {
        // Partial backward-compatibility
        // Converts the regions entry in the regionProperties dictionary into
        // the regionSolvers list

        typeIOobject<IOdictionary> regionPropertiesHeader
        (
            IOobject
            (
                "regionProperties",
                runTime.time().constant(),
                runTime.db(),
                IOobject::MUST_READ
            )
        );

        if (regionPropertiesHeader.headerOk())
        {
            HashTable<wordList> regions
            (
                IOdictionary(regionPropertiesHeader).lookup("regions")
            );

            if (regions.found("solid"))
            {
                const wordList& fluidRegions = regions["solid"];
                forAll(fluidRegions, i)
                {
                    append(Pair<word>(fluidRegions[i], "solid"));
                }
            }

            if (regions.found("fluid"))
            {
                const wordList& fluidRegions = regions["fluid"];
                forAll(fluidRegions, i)
                {
                    append(Pair<word>(fluidRegions[i], "fluid"));
                }
            }
        }
        else
        {
            FatalIOErrorInFunction(runTime.controlDict())
                << "regionSolvers list missing from "
                << runTime.controlDict().name()
                << exit(FatalIOError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolvers::~regionSolvers()
{}


// ************************************************************************* //
