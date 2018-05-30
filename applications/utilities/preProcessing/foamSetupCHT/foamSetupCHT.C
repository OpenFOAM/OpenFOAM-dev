/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

    The utility reads constant/materialProperties to create a
    regionProperties file and to create region directories containing
    required files within the 0, system and constant directories. The
    materialProperties file contains mesh region names with an
    associated physical region type and a material:

    bottomAir
    {
        type            fluid;
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
      + air: radiationProperties, thermophysicalProperties, turbulenceProperties
      + aluminium: radiationProperties, thermophysicalProperties
      + ...

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"

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

    IOdictionary regionProperties
    (
        IOobject
        (
            "regionProperties",
            runTime.constant(),
            runTime
        ),
        dictionary::null
    );

    const fileName currentDir(runTime.path());
    const fileName materialsDir(currentDir/"templates"/"materials");
    const fileName systemDir(currentDir/"templates"/"system");
    const fileName constantDir(currentDir/"templates"/"constant");
    const fileName timeDir(currentDir/"templates"/"0");

    HashTable<wordList> regionInfo;

    forAllConstIter(dictionary, materialProperties, regionIter)
    {
        // Read region subdict name, then its type and material entries
        const word& regionName = regionIter().keyword();
        const dictionary& regionDict = regionIter().dict();
        const word regionType(regionDict.lookup("type"));
        const word regionMaterial(regionDict.lookup("material"));

        Info<< "\nRegion " << regionName << ":\n"
            << "\tCreating 0/" << regionName
            << " directory" << endl;

        // Create fluid and solid lists for regionProperties, e.g.
        // fluid ( bottomAir ... )
        HashTable<wordList>::iterator iter = regionInfo.find(regionType);

        if (iter == regionInfo.end())
        {
            // Add fluid or solid list if does not exist
            regionInfo.insert(regionType, wordList(1, regionName));
        }
        else
        {
            wordList& regionNames = iter();
            if (findIndex(regionNames, regionName) == -1)
            {
                // Append region name to fluid/solid list
                regionNames.append(regionName);
            }
        }

        // 0/<region>: from fluid/solid template
        const fileName sourceDir(timeDir/regionType);
        if (isDir(sourceDir))
        {
            cpFiles(sourceDir, currentDir/"0"/regionName);
        }
        else
        {
            // This needs checking ...
            FatalIOErrorIn(args.executable().c_str(), materialProperties)
                << "Cannot find region type file "
                << sourceDir << exit(FatalIOError);
        }

        Info<< "\tCreating constant/" << regionName
            << " directory with " << regionMaterial
            << " material" << endl;
        cpFiles(constantDir/regionType, currentDir/"constant"/regionName);
        cpFiles(materialsDir/regionMaterial, currentDir/"constant"/regionName);

        // system/<region>: from fluid or solid templ
        Info<< "\tCreating system/" << regionName
            << " directory" << endl;
        cpFiles(systemDir/regionType, currentDir/"system"/regionName);
    }

    regionProperties.add("regions", regionInfo);

    Info<< "\nWriting region properties\n" << endl;

    regionProperties.regIOobject::write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
