/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
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

    dictionaryEntry regionSolvers
    (
        "regionSolvers",
        dictionary::null,
        dictionary::null
    );

    forAllConstIter(dictionary, materialProperties, regionIter)
    {
        // Read region subdict name, then its type and material entries
        const word& regionName = regionIter().keyword();
        const dictionary& regionDict = regionIter().dict();
        const word regionSolver(regionDict.lookup("solver"));
        const word regionMaterial(regionDict.lookup("material"));

        Info<< "\nRegion " << regionName << ":\n"
            << "\tCreating 0/" << regionName
            << " directory" << endl;

        regionSolvers.add(regionName, regionSolver);

        // 0/<region>: from fluid/solid template
        const fileName sourceDir(timeDir/regionSolver);
        if (isDir(sourceDir))
        {
            cpFiles(sourceDir, currentDir/"0"/regionName);
        }
        else
        {
            FatalIOErrorIn(args.executable().c_str(), materialProperties)
                << "Cannot find region type file "
                << sourceDir << exit(FatalIOError);
        }

        const fileName matDir(materialsDir/regionMaterial);
        if (isDir(matDir))
        {
            Info<< "\tCreating constant/" << regionName
                << " directory with " << regionMaterial
                << " material" << endl;
            cpFiles(constantDir/regionSolver, currentDir/"constant"/regionName);
            cpFiles(matDir, currentDir/"constant"/regionName);

            // system/<region>: from fluid or solid template
            Info<< "\tCreating system/" << regionName
                << " directory" << endl;
            cpFiles(systemDir/regionSolver, currentDir/"system"/regionName);
        }
        else
        {
            FatalIOErrorIn(args.executable().c_str(), materialProperties)
                << "Cannot find region material folder "
                << regionMaterial << exit(FatalIOError);
        }
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
