/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "writeMeshObject.H"
#include "fieldDictionary.H"

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

    if (io.headerOk())
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
            runTime.writeCompression()
        );
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

#   include "addRegionOption.H"
#   include "setRootCase.H"

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


#   include "createTime.H"


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
        writeMeshObject<cellCompactIOList>("cells", meshDir, runTime);
        writeMeshObject<labelIOList>("owner", meshDir, runTime);
        writeMeshObject<labelIOList>("neighbour", meshDir, runTime);
        writeMeshObject<faceCompactIOList>("faces", meshDir, runTime);
        writeMeshObject<pointIOField>("points", meshDir, runTime);
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
            )
            {
                Info<< "        Reading " << headerClassName
                    << " : " << iter()->name() << endl;

                fieldDictionary fDict
                (
                    *iter(),
                    headerClassName
                );

                Info<< "        Writing " << iter()->name() << endl;
                fDict.regIOobject::write();
            }
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
