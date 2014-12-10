/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    foamToEnsightParts

Description
    Translates OpenFOAM data to Ensight format.
    An Ensight part is created for each cellZone and patch.

Usage
    - foamToEnsightParts [OPTION] \n
    Translates OpenFOAM data to Ensight format

    \param -ascii \n
    Write Ensight data in ASCII format instead of "C Binary"

    \parm -name \<subdir\>\n
    define sub-directory name to use for Ensight data (default: "Ensight")

    \param -noZero \n
    Exclude the often incomplete initial conditions.

    \param -index \<start\>\n
    Ignore the time index contained in the time file and use a
    simple indexing when creating the \c Ensight/data/######## files.

    \param -noMesh \n
    Suppress writing the geometry. Can be useful for converting partial
    results for a static geometry.

    \param -width \<n\>\n
    width of Ensight data subdir

Note
    - no parallel data.
    - writes to \a Ensight directory to avoid collisions with foamToEnsight.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "volFields.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "IOobjectList.H"
#include "scalarIOField.H"
#include "tensorIOField.H"

#include "ensightParts.H"
#include "ensightOutputFunctions.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    // enable -constant
    // probably don't need -zeroTime though, since the fields are vetted
    // afterwards anyhow
    timeSelector::addOptions(true, false);
    argList::noParallel();
    argList::addBoolOption
    (
        "ascii",
        "write in ASCII format instead of 'C Binary'"
    );
    argList::addOption
    (
        "index",
        "start",
        "ignore the time index contained in the uniform/time file "
        "and use simple indexing when creating the files"
    );
    argList::addBoolOption
    (
        "noMesh",
        "suppress writing the geometry. "
        "Can be useful for converting partial results for a static geometry"
    );
    argList::addOption
    (
        "name",
        "subdir",
        "define sub-directory name to use for Ensight data "
        "(default: \"Ensight\")"
    );
    argList::addOption
    (
        "width",
        "n",
        "width of Ensight data subdir"
    );

    // the volume field types that we handle
    wordHashSet volFieldTypes;
    volFieldTypes.insert(volScalarField::typeName);
    volFieldTypes.insert(volVectorField::typeName);
    volFieldTypes.insert(volSphericalTensorField::typeName);
    volFieldTypes.insert(volSymmTensorField::typeName);
    volFieldTypes.insert(volTensorField::typeName);

    // the lagrangian field types that we handle
    wordHashSet cloudFieldTypes;
    cloudFieldTypes.insert(scalarIOField::typeName);
    cloudFieldTypes.insert(vectorIOField::typeName);
    cloudFieldTypes.insert(tensorIOField::typeName);

    const char* geometryName = "geometry";

    #include "setRootCase.H"
    #include "createTime.H"

    // get times list
    instantList timeDirs = timeSelector::select0(runTime, args);

    // default to binary output, unless otherwise specified
    IOstream::streamFormat format = IOstream::BINARY;
    if (args.optionFound("ascii"))
    {
        format = IOstream::ASCII;
    }

    // control for renumbering iterations
    label indexingNumber = 0;
    bool optIndex = args.optionReadIfPresent("index", indexingNumber);

    // always write the geometry, unless the -noMesh option is specified
    bool optNoMesh = args.optionFound("noMesh");

    // adjust output width
    if (args.optionFound("width"))
    {
        ensightFile::subDirWidth(args.optionRead<label>("width"));
    }

    // define sub-directory name to use for Ensight data
    fileName ensightDir = "Ensight";
    args.optionReadIfPresent("name", ensightDir);

    if (!ensightDir.isAbsolute())
    {
        ensightDir = args.rootPath()/args.globalCaseName()/ensightDir;
    }

    fileName dataDir = ensightDir/"data";
    fileName caseFileName = "Ensight.case";
    fileName dataMask = fileName("data")/ensightFile::mask();

    // Ensight and Ensight/data directories must exist
    // do not remove old data - we might wish to convert new results
    // or a particular time interval
    if (isDir(ensightDir))
    {
        Info<<"Warning: re-using existing directory" << nl
            << "    " << ensightDir << endl;
    }
    mkDir(ensightDir);
    mkDir(dataDir);

    #include "createNamedMesh.H"

    // Mesh instance (region0 gets filtered out)
    fileName regionPrefix;

    if (regionName != polyMesh::defaultRegion)
    {
        regionPrefix = regionName;
    }

    // Construct the list of ensight parts for the entire mesh
    ensightParts partsList(mesh);

    // write summary information
    {
        OFstream partsInfoFile(ensightDir/"partsInfo");

        partsInfoFile
            << "// summary of ensight parts" << nl << nl;
        partsList.writeSummary(partsInfoFile);
    }

    #include "checkHasMovingMesh.H"
    #include "findFields.H"

    if (hasMovingMesh && optNoMesh)
    {
        Info<< "mesh is moving: ignoring '-noMesh' option" << endl;
        optNoMesh = false;
    }


    // map times used
    Map<scalar>  timeIndices;

    // Track the time indices used by the volume fields
    DynamicList<label> fieldTimesUsed;

    // Track the time indices used by each cloud
    HashTable<DynamicList<label> > cloudTimesUsed;

    // Create a new DynamicList for each cloud
    forAllConstIter(HashTable<HashTable<word> >, cloudFields, cloudIter)
    {
        cloudTimesUsed.insert(cloudIter.key(), DynamicList<label>());
    }


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        #include "getTimeIndex.H"

        // remember the time index
        fieldTimesUsed.append(timeIndex);

        // the data/ITER subdirectory must exist
        fileName subDir = ensightFile::subDir(timeIndex);
        mkDir(dataDir/subDir);

        // place a timestamp in the directory for future reference
        {
            OFstream timeStamp(dataDir/subDir/"time");
            timeStamp
                << "#   timestep time" << nl
                << subDir.c_str() << " " << runTime.timeName() << nl;
        }

        #include "moveMesh.H"

        if (timeI == 0 || mesh.moving())
        {
            if (mesh.moving())
            {
                partsList.recalculate(mesh);
            }

            if (!optNoMesh)
            {
                fileName geomDir;
                if (hasMovingMesh)
                {
                    geomDir = dataDir/subDir;
                }

                ensightGeoFile geoFile(ensightDir/geomDir/geometryName, format);
                partsList.writeGeometry(geoFile);
                Info<< nl;
            }
        }

        Info<< "write volume field (" << flush;

        forAllConstIter(HashTable<word>, volumeFields, fieldIter)
        {
            const word& fieldName = fieldIter.key();
            const word& fieldType = fieldIter();

            IOobject fieldObject
            (
                fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            );

            if (fieldType == volScalarField::typeName)
            {
                ensightVolField<scalar>
                (
                    partsList,
                    fieldObject,
                    mesh,
                    dataDir,
                    subDir,
                    format
                );

            }
            else if (fieldType == volVectorField::typeName)
            {
                ensightVolField<vector>
                (
                    partsList,
                    fieldObject,
                    mesh,
                    dataDir,
                    subDir,
                    format
                );

            }
            else if (fieldType == volSphericalTensorField::typeName)
            {
                ensightVolField<sphericalTensor>
                (
                    partsList,
                    fieldObject,
                    mesh,
                    dataDir,
                    subDir,
                    format
                );

            }
            else if (fieldType == volSymmTensorField::typeName)
            {
                ensightVolField<symmTensor>
                (
                    partsList,
                    fieldObject,
                    mesh,
                    dataDir,
                    subDir,
                    format
                );
            }
            else if (fieldType == volTensorField::typeName)
            {
                ensightVolField<tensor>
                (
                    partsList,
                    fieldObject,
                    mesh,
                    dataDir,
                    subDir,
                    format
                );
            }
        }
        Info<< " )" << endl;

        // check for clouds
        forAllConstIter(HashTable<HashTable<word> >, cloudFields, cloudIter)
        {
            const word& cloudName = cloudIter.key();

            if
            (
                !isDir
                (
                    runTime.timePath()/regionPrefix/
                    cloud::prefix/cloudName
                )
            )
            {
                continue;
            }

            IOobjectList cloudObjs
            (
                mesh,
                runTime.timeName(),
                cloud::prefix/cloudName
            );

            // check that the positions field is present for this time
            IOobject* positionPtr = cloudObjs.lookup(word("positions"));
            if (positionPtr != NULL)
            {
                ensightParticlePositions
                (
                    mesh,
                    dataDir,
                    subDir,
                    cloudName,
                    format
                );
            }
            else
            {
                continue;
            }

            Info<< "write " << cloudName << " (" << flush;

            forAllConstIter(HashTable<word>, cloudIter(), fieldIter)
            {
                const word& fieldName = fieldIter.key();
                const word& fieldType = fieldIter();

                IOobject *fieldObject = cloudObjs.lookup(fieldName);

                if (!fieldObject)
                {
                    Info<< "missing "
                        << runTime.timeName()/cloud::prefix/cloudName
                        / fieldName
                        << endl;
                    continue;
                }

                if (fieldType == scalarIOField::typeName)
                {
                    ensightLagrangianField<scalar>
                    (
                        *fieldObject,
                        dataDir,
                        subDir,
                        cloudName,
                        format
                    );

                }
                else if (fieldType == vectorIOField::typeName)
                {
                    ensightLagrangianField<vector>
                    (
                        *fieldObject,
                        dataDir,
                        subDir,
                        cloudName,
                        format
                    );

                }
                else if (fieldType == tensorIOField::typeName)
                {
                    ensightLagrangianField<tensor>
                    (
                        *fieldObject,
                        dataDir,
                        subDir,
                        cloudName,
                        format
                    );

                }
            }

            Info<< " )" << endl;

            // remember the time index
            cloudTimesUsed[cloudName].append(timeIndex);
        }
    }

    #include "ensightOutputCase.H"

    Info<< "\nEnd\n"<< endl;

    return 0;
}


// ************************************************************************* //
