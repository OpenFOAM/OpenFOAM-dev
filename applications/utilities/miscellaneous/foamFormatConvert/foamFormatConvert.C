/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    foamFormatConvert

Description
    Converts all IOobjects associated with a case into the format specified
    in the controlDict.

    Mainly used to convert binary mesh/field files to ASCII.

    Problem: any zero-size List written binary gets written as '0'. When
    reading the file as a dictionary this is interpreted as a label. This
    is (usually) not a problem when doing patch fields since these get the
    'uniform', 'nonuniform' prefix. However zone contents are labelLists
    not labelFields and these go wrong. For now hacked a solution where
    we detect the keywords in zones and redo the dictionary entries
    to be labelLists.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "cellIOList.H"
#include "IOobjectList.H"
#include "IOPtrList.H"
#include "cloud.H"
#include "labelIOField.H"
#include "scalarIOField.H"
#include "sphericalTensorIOField.H"
#include "symmTensorIOField.H"
#include "tensorIOField.H"
#include "labelFieldIOField.H"
#include "vectorFieldIOField.H"
#include "Cloud.H"
#include "passiveParticle.H"
#include "fieldDictionary.H"

#include "writeMeshObject.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<entry>, 0);
}


// Hack to do zones which have Lists in them. See above.
bool writeZones(const word& name, const fileName& meshDir, Time& runTime)
{
    IOobject io
    (
        name,
        runTime.timeName(),
        meshDir,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    bool writeOk = false;

    if (io.typeHeaderOk<meshCellZones>(false))
    {
        Info<< "        Reading " << io.headerClassName()
            << " : " << name << endl;

        // Switch off type checking (for reading e.g. faceZones as
        // generic list of dictionaries).
        const word oldTypeName = IOPtrList<entry>::typeName;
        const_cast<word&>(IOPtrList<entry>::typeName) = word::null;

        IOPtrList<entry> meshObject(io);

        forAll(meshObject, i)
        {
            if (meshObject[i].isDict())
            {
                dictionary& d = meshObject[i].dict();

                if (d.found("faceLabels"))
                {
                    d.set("faceLabels", labelList(d.lookup("faceLabels")));
                }

                if (d.found("flipMap"))
                {
                    d.set("flipMap", boolList(d.lookup("flipMap")));
                }

                if (d.found("cellLabels"))
                {
                    d.set("cellLabels", labelList(d.lookup("cellLabels")));
                }

                if (d.found("pointLabels"))
                {
                    d.set("pointLabels", labelList(d.lookup("pointLabels")));
                }
            }
        }

        const_cast<word&>(IOPtrList<entry>::typeName) = oldTypeName;
        // Fake type back to what was in field
        const_cast<word&>(meshObject.type()) = io.headerClassName();

        Info<< "        Writing " << name << endl;

        // Force writing as ascii
        writeOk = meshObject.regIOobject::writeObject
        (
            IOstream::ASCII,
            IOstream::currentVersion,
            runTime.writeCompression(),
            true
        );
    }

    return writeOk;
}


// Reduction for non-empty strings
class uniqueEqOp
{
    public:
    void operator()(stringList& x, const stringList& y) const
    {
        stringList newX(x.size()+y.size());
        label n = 0;
        forAll(x, i)
        {
            if (!x[i].empty())
            {
                newX[n++] = x[i];
            }
        }
        forAll(y, i)
        {
            if (!y[i].empty() && findIndex(x, y[i]) == -1)
            {
                newX[n++] = y[i];
            }
        }
        newX.setSize(n);
        x.transfer(newX);
    }
};


template<class T>
bool writeOptionalMeshObject
(
    const word& name,
    const fileName& meshDir,
    Time& runTime,
    const bool write
)
{
    IOobject io
    (
        name,
        runTime.timeName(),
        meshDir,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    bool writeOk = false;

    bool haveFile = io.typeHeaderOk<IOField<label>>(false);

    // Make sure all know if there is a valid class name
    stringList classNames(1, io.headerClassName());
    combineReduce(classNames, uniqueEqOp());

    // Check for correct type
    if (classNames[0] == T::typeName)
    {
        Info<< "        Reading " << classNames[0]
            << " : " << name << endl;
        T meshObject(io, write && haveFile);

        Info<< "        Writing " << name << endl;
        writeOk = meshObject.regIOobject::write(write && haveFile);
    }

    return writeOk;
}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::addBoolOption
    (
        "noConstant",
        "exclude the 'constant/' dir in the times list"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"

    // enable noConstant by switching
    if (!args.optionFound("noConstant"))
    {
        args.setOption("constant", "");
    }
    else
    {
        args.unsetOption("constant");
        Info<< "Excluding the constant directory." << nl << endl;
    }


    #include "createTime.H"
    // Optional mesh (used to read Clouds)
    autoPtr<polyMesh> meshPtr;


    // Make sure we do not use the master-only reading since we read
    // fields (different per processor) as dictionaries.
    regIOobject::fileModificationChecking = regIOobject::timeStamp;


    fileName meshDir = polyMesh::meshSubDir;
    fileName regionPrefix = "";
    word regionName = polyMesh::defaultRegion;
    if (args.optionReadIfPresent("region", regionName))
    {
        Info<< "Using region " << regionName << nl << endl;
        regionPrefix = regionName;
        meshDir = regionName/polyMesh::meshSubDir;
    }

    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        // Convert all the standard mesh files
        writeMeshObject<cellCompactIOList, cellIOList>
        (
            "cells",
            meshDir,
            runTime
        );
        writeMeshObject<labelIOList>("owner", meshDir, runTime);
        writeMeshObject<labelIOList>("neighbour", meshDir, runTime);
        writeMeshObject<faceCompactIOList, faceIOList>
        (
            "faces",
            meshDir,
            runTime
        );
        writeMeshObject<pointIOField>("points", meshDir, runTime);
        // Write boundary in ascii. This is only needed for fileHandler to
        // kick in. Should not give problems since always writing ascii.
        writeZones("boundary", meshDir, runTime);
        writeMeshObject<labelIOList>("pointProcAddressing", meshDir, runTime);
        writeMeshObject<labelIOList>("faceProcAddressing", meshDir, runTime);
        writeMeshObject<labelIOList>("cellProcAddressing", meshDir, runTime);
        writeMeshObject<labelIOList>
        (
            "boundaryProcAddressing",
            meshDir,
            runTime
        );

        // foamyHexMesh vertices
        writeMeshObject<pointIOField>
        (
            "internalDelaunayVertices",
            regionPrefix,
            runTime
        );

        if (runTime.writeFormat() == IOstream::ASCII)
        {
            // Only do zones when converting from binary to ascii
            // The other way gives problems since working on dictionary level.
            writeZones("cellZones", meshDir, runTime);
            writeZones("faceZones", meshDir, runTime);
            writeZones("pointZones", meshDir, runTime);
        }

        // Get list of objects from the database
        IOobjectList objects(runTime, runTime.timeName(), regionPrefix);

        forAllConstIter(IOobjectList, objects, iter)
        {
            const word& headerClassName = iter()->headerClassName();

            if
            (
                headerClassName == volScalarField::typeName
             || headerClassName == volVectorField::typeName
             || headerClassName == volSphericalTensorField::typeName
             || headerClassName == volSymmTensorField::typeName
             || headerClassName == volTensorField::typeName

             || headerClassName == surfaceScalarField::typeName
             || headerClassName == surfaceVectorField::typeName
             || headerClassName == surfaceSphericalTensorField::typeName
             || headerClassName == surfaceSymmTensorField::typeName
             || headerClassName == surfaceTensorField::typeName

             || headerClassName == pointScalarField::typeName
             || headerClassName == pointVectorField::typeName
             || headerClassName == pointSphericalTensorField::typeName
             || headerClassName == pointSymmTensorField::typeName
             || headerClassName == pointTensorField::typeName

             || headerClassName == volScalarField::Internal::typeName
             || headerClassName == volVectorField::Internal::typeName
             || headerClassName == volSphericalTensorField::Internal::typeName
             || headerClassName == volSymmTensorField::Internal::typeName
             || headerClassName == volTensorField::Internal::typeName
            )
            {
                Info<< "        Reading " << headerClassName
                    << " : " << iter()->name() << endl;

                fieldDictionary fDict(*iter(), headerClassName);

                Info<< "        Writing " << iter()->name() << endl;
                fDict.regIOobject::write();
            }
        }



        // Check for lagrangian
        const fileName lagrangianDir
        (
            fileHandler().filePath
            (
                runTime.timePath()
              / regionPrefix
              / cloud::prefix
            )
        );
        stringList lagrangianDirs
        (
            lagrangianDir == fileName::null ? 0 : 1,
            lagrangianDir
        );

        combineReduce(lagrangianDirs, uniqueEqOp());

        if (!lagrangianDirs.empty())
        {
            if (meshPtr.valid())
            {
                meshPtr().readUpdate();
            }
            else
            {
                Info<< "        Create polyMesh for time = "
                    << runTime.timeName() << endl;

                meshPtr.reset
                (
                    new polyMesh
                    (
                        IOobject
                        (
                            regionName,
                            runTime.timeName(),
                            runTime,
                            Foam::IOobject::MUST_READ
                        )
                    )
                );
            }

            stringList cloudDirs
            (
                fileHandler().readDir
                (
                    lagrangianDirs[0],
                    fileType::directory
                )
            );

            combineReduce(cloudDirs, uniqueEqOp());

            forAll(cloudDirs, i)
            {
                fileName dir(cloud::prefix/cloudDirs[i]);

                Cloud<passiveParticle> parcels(meshPtr(), cloudDirs[i], false);

                parcels.writeObject
                (
                    runTime.writeFormat(),
                    IOstream::currentVersion,
                    runTime.writeCompression(),
                    parcels.size()
                );


                // Do local scan for valid cloud objects
                IOobjectList sprayObjs(runTime, runTime.timeName(), dir);

                // Combine with all other cloud objects
                stringList sprayFields(sprayObjs.sortedToc());
                combineReduce(sprayFields, uniqueEqOp());

                forAll(sprayFields, fieldi)
                {
                    const word& name = sprayFields[fieldi];

                    // Note: try the various field types. Make sure to
                    //       exit once successful conversion to avoid re-read
                    //       converted file.

                    if
                    (
                        name == "positions"
                     || name == "origProcId"
                     || name == "origId"
                    )
                    {
                        continue;
                    }

                    bool writeOk = writeOptionalMeshObject<labelIOField>
                    (
                        name,
                        dir,
                        runTime,
                        parcels.size() > 0
                    );
                    if (writeOk) continue;

                    writeOk = writeOptionalMeshObject<scalarIOField>
                    (
                        name,
                        dir,
                        runTime,
                        parcels.size() > 0
                    );
                    if (writeOk) continue;

                    writeOk = writeOptionalMeshObject<vectorIOField>
                    (
                        name,
                        dir,
                        runTime,
                        parcels.size() > 0
                    );
                    if (writeOk) continue;

                    writeOk = writeOptionalMeshObject<sphericalTensorIOField>
                    (
                        name,
                        dir,
                        runTime,
                        parcels.size() > 0
                    );
                    if (writeOk) continue;

                    writeOk = writeOptionalMeshObject<symmTensorIOField>
                    (
                        name,
                        dir,
                        runTime,
                        parcels.size() > 0
                    );
                    if (writeOk) continue;

                    writeOk = writeOptionalMeshObject<tensorIOField>
                    (
                        name,
                        dir,
                        runTime,
                        parcels.size() > 0
                    );
                    if (writeOk) continue;

                    writeOk = writeOptionalMeshObject<labelFieldIOField>
                    (
                        name,
                        dir,
                        runTime,
                        parcels.size() > 0
                    );
                    if (writeOk) continue;

                    writeOk = writeOptionalMeshObject<vectorFieldIOField>
                    (
                        name,
                        dir,
                        runTime,
                        parcels.size() > 0
                    );

                    if (!writeOk)
                    {
                        Info<< "        Failed converting " << name << endl;
                    }
                }
            }
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
