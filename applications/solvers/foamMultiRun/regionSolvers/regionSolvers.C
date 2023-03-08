/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
#include "solver.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolvers::regionSolvers(const Time& runTime)
{
    List<Pair<word>> regionSolverNames;

    if (runTime.controlDict().found("regionSolvers"))
    {
        const dictionary& regionSolversDict =
            runTime.controlDict().subDict("regionSolvers");

        forAllConstIter(dictionary, regionSolversDict, iter)
        {
            const word regionName(iter().keyword());
            const word solverName(iter().stream());

            regionSolverNames.append(Pair<word>(regionName, solverName));
        }
    }
    else
    {
        // Partial backward-compatibility
        // Converts the regions entry in the regionProperties dictionary into
        // the regionSolvers list
        // Only supports fluid and solid regions

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
                    regionSolverNames.append
                    (
                        Pair<word>(fluidRegions[i], "solid")
                    );
                }
            }

            if (regions.found("fluid"))
            {
                const wordList& fluidRegions = regions["fluid"];
                forAll(fluidRegions, i)
                {
                    regionSolverNames.append
                    (
                        Pair<word>(fluidRegions[i], "fluid")
                    );
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

    regions_.setSize(regionSolverNames.size());
    solvers_.setSize(regionSolverNames.size());
    prefixes_.setSize(regionSolverNames.size());

    string::size_type nRegionNameChars = 0;

    forAll(regionSolverNames, i)
    {
        const word& regionName = regionSolverNames[i].first();
        const word& solverName = regionSolverNames[i].second();

        // Load the solver library
        solver::load(solverName);

        regions_.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    regionName,
                    runTime.name(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );

        solvers_.set(i, solver::New(solverName, regions_[i]));

        prefixes_[i] = regionName;
        nRegionNameChars = max(nRegionNameChars, regionName.size());
    }

    nRegionNameChars++;

    prefix0_.append(nRegionNameChars, ' ');

    forAll(regionSolverNames, i)
    {
        prefixes_[i].append(nRegionNameChars - prefixes_[i].size(), ' ');
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolvers::~regionSolvers()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionSolvers::setGlobalPrefix() const
{
    Sout.prefix() = prefix0_;
}


void Foam::regionSolvers::setPrefix(const label i) const
{
    Sout.prefix() = prefixes_[i];
}


void Foam::regionSolvers::resetPrefix() const
{
    Sout.prefix() = string::null;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::solver& Foam::regionSolvers::operator[](const label i)
{
    setPrefix(i);
    return solvers_[i];
}


// ************************************************************************* //
