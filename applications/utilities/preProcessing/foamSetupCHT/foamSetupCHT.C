/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2025 OpenFOAM Foundation
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

Application
    foamSetupCHT

Description
    This utility sets up a multi-region case using template files for
    material properties, field and system files.

    The utility reads constant/materialProperties to create a regionSolvers list
    and to create region directories containing required files within the 0,
    system and constant directories. The materialProperties file contains mesh
    region names with an associated solver and a material:

    bottomAir
    {
        solver          fluid;
        material        air;
    }

    The case must contain a directory called templates, with e.g. the
    following directories and files:
    + 0
      + fluid: p, p_rgh, U, T, k, omega, epsilon, nut, alphat
      + solid: p, T
    + system
      + fluid: fvSchemes, fvSolution, decomposeParDict
      + solid: fvSchemes, fvSolution, decomposeParDict
    + constant
      + fluid: g
      + solid
    + materials
      + air: radiationProperties, thermophysicalProperties, momentumTransport
      + aluminium: radiationProperties, thermophysicalProperties
      + ...

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "dictionaryEntry.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "Reading material properties\n" << endl;

    IOdictionary materialProperties
    (
        IOobject
        (
            "materialProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const fileName currentDir(runTime.path());
    const fileName materialsDir(currentDir/"templates"/"materials");
    const fileName systemDir(currentDir/"templates"/"system");
    const fileName constantDir(currentDir/"templates"/"constant");
    const fileName timeDir(currentDir/"templates"/"0");

    // Set of valid region names, corresponding to polyMesh sub-directories
    wordHashSet regions
    (
        wordList(readDir(currentDir/"constant", fileType::directory))
    );

    wordList regionDirs(regions.toc());

    forAll(regionDirs, i)
    {
        if(!isDir(currentDir/"constant"/regionDirs[i]/"polyMesh"))
        {
            regions.erase(regionDirs[i]);
        }
    }

    // Set of valid solver names, listed by:
    /*
        grep -rwl thermo $FOAM_MODULES | \
            sed "s@$FOAM_MODULES/\([^/]*\).*@\1@g" | sort -u
    */
    const wordHashSet solvers
    {
        "compressibleMultiphaseVoF",
        "compressibleVoF",
        "film",
        "fluid",
        "isothermalFilm",
        "isothermalFluid",
        "multicomponentFluid",
        "multiphaseEuler",
        "shockFluid",
        "solid",
        "solidDisplacement",
        "XiFluid"
    };

    // Set of valid material names
    wordHashSet materials
    (
        wordList(readDir(materialsDir, fileType::directory))
    );

    wordList regionList;
    wordList solverList;
    wordList materialList;

    forAllConstIter(dictionary, materialProperties, regionIter)
    {
        regionList.append(regionIter().keyword());
        const dictionary& regionDict = regionIter().dict();
        solverList.append(regionDict.lookup("solver"));
        materialList.append(regionDict.lookup("material"));

        if (!regions.found(regionList.last()))
        {
            FatalIOErrorIn(args.executable().c_str(), materialProperties)
                << "Cannot find region polyMesh directory "
                << "constant"/regionList.last()/"polyMesh" << nl << nl
                << "Existing regions are "
                << regions
                << exit(FatalIOError);
        }

        if (!solvers.found(solverList.last()))
        {
            FatalIOErrorIn(args.executable().c_str(), materialProperties)
                << "Cannot find solver "
                << solverList.last() << nl << nl
                << "Available solvers are "
                << solvers.sortedToc()
                << exit(FatalIOError);
        }

        if (!materials.found(materialList.last()))
        {
            FatalIOErrorIn(args.executable().c_str(), materialProperties)
                << "Cannot find region material directory "
                << "templates/materials"/materialList.last() << nl << nl
                << "Available materials are "
                << materials
                << exit(FatalIOError);
        }
    }

    dictionaryEntry regionSolvers
    (
        "regionSolvers",
        dictionary::null,
        dictionary::null
    );

    forAll(regionList, i)
    {
        const word region(regionList[i]);
        const word material(materialList[i]);
        const word solver(solverList[i]);

        Info<< "\nRegion " << region << ":\n"
            << "\tCreating 0/" << region << " directory" << endl;
        cpFiles(timeDir/solver, currentDir/"0"/region);

        Info<< "\tCreating constant/" << region
            << " directory with " << material << " material" << endl;
        cpFiles(constantDir/solver, currentDir/"constant"/region);
        cpFiles(materialsDir/material, currentDir/"constant"/region);

        Info<< "\tCreating system/" << region << " directory" << endl;
        cpFiles(systemDir/solver, currentDir/"system"/region);

        regionSolvers.add(region, solver);
    }

    Info<< "\nWriting regionSolvers\n" << endl;
    {
        OFstream regionSolversFile(currentDir/runTime.system()/"regionSolvers");
        regionSolvers.write(regionSolversFile);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
