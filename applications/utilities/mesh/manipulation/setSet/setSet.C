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
    setSet

Description
    Manipulate a cell/face/point/ set or zone interactively.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "globalMeshData.H"
#include "IStringStream.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "topoSetSource.H"
#include "OFstream.H"
#include "IFstream.H"
#include "demandDrivenData.H"
#include "writePatch.H"
#include "writePointSet.H"
#include "IOobjectList.H"
#include "cellZoneSet.H"
#include "faceZoneSet.H"
#include "pointZoneSet.H"
#include "timeSelector.H"

#include <stdio.h>


#ifdef HAS_READLINE
# include <readline/readline.h>
# include <readline/history.h>
#endif

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef HAS_READLINE
static const char* historyFile = ".setSet";
#endif


// Write set to VTK readable files
void writeVTK
(
    const polyMesh& mesh,
    const topoSet& currentSet,
    const fileName& vtkName
)
{
    if (isA<faceSet>(currentSet))
    {
        // Faces of set with OpenFOAM faceID as value

        faceList setFaces(currentSet.size());
        labelList faceValues(currentSet.size());
        label setFaceI = 0;

        forAllConstIter(topoSet, currentSet, iter)
        {
            setFaces[setFaceI] = mesh.faces()[iter.key()];
            faceValues[setFaceI] = iter.key();
            setFaceI++;
        }

        primitiveFacePatch fp(setFaces, mesh.points());

        writePatch
        (
            true,
            currentSet.name(),
            fp,
            "faceID",
            faceValues,
            mesh.time().path()/vtkName
        );
    }
    else if (isA<cellSet>(currentSet))
    {
        // External faces of cellset with OpenFOAM cellID as value

        Map<label> cellFaces(currentSet.size());

        forAllConstIter(cellSet, currentSet, iter)
        {
            label cellI = iter.key();

            const cell& cFaces = mesh.cells()[cellI];

            forAll(cFaces, i)
            {
                label faceI = cFaces[i];

                if (mesh.isInternalFace(faceI))
                {
                    label otherCellI = mesh.faceOwner()[faceI];

                    if (otherCellI == cellI)
                    {
                        otherCellI = mesh.faceNeighbour()[faceI];
                    }

                    if (!currentSet.found(otherCellI))
                    {
                        cellFaces.insert(faceI, cellI);
                    }
                }
                else
                {
                    cellFaces.insert(faceI, cellI);
                }
            }
        }

        faceList setFaces(cellFaces.size());
        labelList faceValues(cellFaces.size());
        label setFaceI = 0;

        forAllConstIter(Map<label>, cellFaces, iter)
        {
            setFaces[setFaceI] = mesh.faces()[iter.key()];
            faceValues[setFaceI] = iter();              // Cell ID
            setFaceI++;
        }

        primitiveFacePatch fp(setFaces, mesh.points());

        writePatch
        (
            true,
            currentSet.name(),
            fp,
            "cellID",
            faceValues,
            mesh.time().path()/vtkName
        );
    }
    else if (isA<pointSet>(currentSet))
    {
        writePointSet
        (
            true,
            mesh,
            currentSet,
            mesh.time().path()/vtkName
        );
    }
    else
    {
        WarningIn
        (
            "void writeVTK"
            "(const polyMesh& mesh, const topoSet& currentSet,"
            "const fileName& vtkName)"
        )   << "Don't know how to handle set of type " << currentSet.type()
            << endl;
    }
}


void printHelp(Ostream& os)
{
    os  << "Please type 'help', 'list', 'quit', 'time ddd'"
        << " or a set command after prompt." << endl
        << "'list' will show all current cell/face/point sets." << endl
        << "'time ddd' will change the current time." << endl
        << endl
        << "A set command should be of the following form" << endl
        << endl
        << "    cellSet|faceSet|pointSet <setName> <action> <source>"
        << endl
        << endl
        << "The <action> is one of" << endl
        << "    list            - prints the contents of the set" << endl
        << "    clear           - clears the set" << endl
        << "    invert          - inverts the set" << endl
        << "    remove          - remove the set" << endl
        << "    new <source>    - sets to set to the source set" << endl
        << "    add <source>    - adds all elements from the source set" << endl
        << "    delete <source> - deletes      ,," << endl
        << "    subset <source> - combines current set with the source set"
        << endl
        << endl
        << "The sources come in various forms. Type a wrong source"
        << " to see all the types available." << endl
        << endl
        << "Example: pick up all cells connected by point or face to patch"
        << " movingWall" << endl
        << endl
        << "Pick up all faces of patch:" << endl
        << "    faceSet f0 new patchToFace movingWall" << endl
        << "Add faces 0,1,2:" << endl
        << "    faceSet f0 add labelToFace (0 1 2)" << endl
        << "Pick up all points used by faces in faceSet f0:" << endl
        << "    pointSet p0 new faceToPoint f0 all" << endl
        << "Pick up cell which has any face in f0:" << endl
        << "    cellSet c0 new faceToCell f0 any" << endl
        << "Add cells which have any point in p0:" << endl
        << "    cellSet c0 add pointToCell p0 any" << endl
        << "List set:" << endl
        << "    cellSet c0 list" << endl
        << endl
        << "Zones can be set using zoneSets from corresponding sets:" << endl
        << "    cellZoneSet c0Zone new setToCellZone c0" << endl
        << "    faceZoneSet f0Zone new setToFaceZone f0" << endl
        << endl
        << "or if orientation is important:" << endl
        << "    faceZoneSet f0Zone new setsToFaceZone f0 c0" << endl
        << endl
        << "ZoneSets can be manipulated using the general actions:" << endl
        << "    list            - prints the contents of the set" << endl
        << "    clear           - clears the set" << endl
        << "    invert          - inverts the set (undefined orientation)"
        << endl
        << "    remove          - remove the set" << endl
        << endl;
}


void printAllSets(const polyMesh& mesh, Ostream& os)
{
    IOobjectList objects
    (
        mesh,
        mesh.time().findInstance
        (
            polyMesh::meshSubDir/"sets",
            word::null,
            IOobject::READ_IF_PRESENT,
            mesh.facesInstance()
        ),
        polyMesh::meshSubDir/"sets"
    );
    IOobjectList cellSets(objects.lookupClass(cellSet::typeName));
    if (cellSets.size())
    {
        os  << "cellSets:" << endl;
        forAllConstIter(IOobjectList, cellSets, iter)
        {
            cellSet set(*iter());
            os  << '\t' << set.name() << "\tsize:" << set.size() << endl;
        }
    }
    IOobjectList faceSets(objects.lookupClass(faceSet::typeName));
    if (faceSets.size())
    {
        os  << "faceSets:" << endl;
        forAllConstIter(IOobjectList, faceSets, iter)
        {
            faceSet set(*iter());
            os  << '\t' << set.name() << "\tsize:" << set.size() << endl;
        }
    }
    IOobjectList pointSets(objects.lookupClass(pointSet::typeName));
    if (pointSets.size())
    {
        os  << "pointSets:" << endl;
        forAllConstIter(IOobjectList, pointSets, iter)
        {
            pointSet set(*iter());
            os  << '\t' << set.name() << "\tsize:" << set.size() << endl;
        }
    }

    const cellZoneMesh& cellZones = mesh.cellZones();
    if (cellZones.size())
    {
        os  << "cellZones:" << endl;
        forAll(cellZones, i)
        {
            const cellZone& zone = cellZones[i];
            os  << '\t' << zone.name() << "\tsize:" << zone.size() << endl;
        }
    }
    const faceZoneMesh& faceZones = mesh.faceZones();
    if (faceZones.size())
    {
        os  << "faceZones:" << endl;
        forAll(faceZones, i)
        {
            const faceZone& zone = faceZones[i];
            os  << '\t' << zone.name() << "\tsize:" << zone.size() << endl;
        }
    }
    const pointZoneMesh& pointZones = mesh.pointZones();
    if (pointZones.size())
    {
        os  << "pointZones:" << endl;
        forAll(pointZones, i)
        {
            const pointZone& zone = pointZones[i];
            os  << '\t' << zone.name() << "\tsize:" << zone.size() << endl;
        }
    }

    os  << endl;
}


template<class ZoneType>
void removeZone
(
    ZoneMesh<ZoneType, polyMesh>& zones,
    const word& setName
)
{
    label zoneID = zones.findZoneID(setName);

    if (zoneID != -1)
    {
        Info<< "Removing zone " << setName << " at index " << zoneID << endl;
        // Shuffle to last position
        labelList oldToNew(zones.size());
        label newI = 0;
        forAll(oldToNew, i)
        {
            if (i != zoneID)
            {
                oldToNew[i] = newI++;
            }
        }
        oldToNew[zoneID] = newI;
        zones.reorder(oldToNew);
        // Remove last element
        zones.setSize(zones.size()-1);
        zones.clearAddressing();
        zones.write();
    }
}


// Physically remove a set
void removeSet
(
    const polyMesh& mesh,
    const word& setType,
    const word& setName
)
{
    // Remove the file
    IOobjectList objects
    (
        mesh,
        mesh.time().findInstance
        (
            polyMesh::meshSubDir/"sets",
            word::null,
            IOobject::READ_IF_PRESENT,
            mesh.facesInstance()
        ),
        polyMesh::meshSubDir/"sets"
    );

    if (objects.found(setName))
    {
        // Remove file
        fileName object = objects[setName]->objectPath();
        Info<< "Removing file " << object << endl;
        rm(object);
    }

    // See if zone
    if (setType == cellZoneSet::typeName)
    {
        removeZone
        (
            const_cast<cellZoneMesh&>(mesh.cellZones()),
            setName
        );
    }
    else if (setType == faceZoneSet::typeName)
    {
        removeZone
        (
            const_cast<faceZoneMesh&>(mesh.faceZones()),
            setName
        );
    }
    else if (setType == pointZoneSet::typeName)
    {
        removeZone
        (
            const_cast<pointZoneMesh&>(mesh.pointZones()),
            setName
        );
    }
}


// Read command and execute. Return true if ok, false otherwise.
bool doCommand
(
    const polyMesh& mesh,
    const word& setType,
    const word& setName,
    const word& actionName,
    const bool writeVTKFile,
    const bool writeCurrentTime,
    const bool noSync,
    Istream& is
)
{
    // Get some size estimate for set.
    const globalMeshData& parData = mesh.globalData();

    label typSize =
        max
        (
            parData.nTotalCells(),
            max
            (
                parData.nTotalFaces(),
                parData.nTotalPoints()
            )
        )
      / (10*Pstream::nProcs());


    bool ok = true;

    // Set to work on
    autoPtr<topoSet> currentSetPtr;

    word sourceType;

    try
    {
        topoSetSource::setAction action =
            topoSetSource::toAction(actionName);


        IOobject::readOption r;

        if (action == topoSetSource::REMOVE)
        {
            removeSet(mesh, setType, setName);
        }
        else if
        (
            (action == topoSetSource::NEW)
         || (action == topoSetSource::CLEAR)
        )
        {
            r = IOobject::NO_READ;
            currentSetPtr = topoSet::New(setType, mesh, setName, typSize);
        }
        else
        {
            r = IOobject::MUST_READ;
            currentSetPtr = topoSet::New(setType, mesh, setName, r);
            topoSet& currentSet = currentSetPtr();
            // Presize it according to current mesh data.
            currentSet.resize(max(currentSet.size(), typSize));
        }

        if (currentSetPtr.valid())
        {
            topoSet& currentSet = currentSetPtr();

            Info<< "    Set:" << currentSet.name()
                << "  Size:" << returnReduce(currentSet.size(), sumOp<label>())
                << "  Action:" << actionName
                << endl;

            switch (action)
            {
                case topoSetSource::CLEAR:
                {
                    // Already handled above by not reading
                    break;
                }
                case topoSetSource::INVERT:
                {
                    currentSet.invert(currentSet.maxSize(mesh));
                    break;
                }
                case topoSetSource::LIST:
                {
                    currentSet.writeDebug(Pout, mesh, 100);
                    Pout<< endl;
                    break;
                }
                case topoSetSource::SUBSET:
                {
                    if (is >> sourceType)
                    {
                        autoPtr<topoSetSource> setSource
                        (
                            topoSetSource::New
                            (
                                sourceType,
                                mesh,
                                is
                            )
                        );

                        // Backup current set.
                        autoPtr<topoSet> oldSet
                        (
                            topoSet::New
                            (
                                setType,
                                mesh,
                                currentSet.name() + "_old2",
                                currentSet
                            )
                        );

                        currentSet.clear();
                        setSource().applyToSet(topoSetSource::NEW, currentSet);

                        // Combine new value of currentSet with old one.
                        currentSet.subset(oldSet);
                    }
                    break;
                }
                default:
                {
                    if (is >> sourceType)
                    {
                        autoPtr<topoSetSource> setSource
                        (
                            topoSetSource::New
                            (
                                sourceType,
                                mesh,
                                is
                            )
                        );

                        setSource().applyToSet(action, currentSet);
                    }
                }
            }


            if (action != topoSetSource::LIST)
            {
                // Set will have been modified.

                // Synchronize for coupled patches.
                if (!noSync) currentSet.sync(mesh);

                // Write
                if (writeVTKFile)
                {
                    mkDir(mesh.time().path()/"VTK"/currentSet.name());

                    fileName vtkName
                    (
                        "VTK"/currentSet.name()/currentSet.name()
                      + "_"
                      + name(mesh.time().timeIndex())
                      + ".vtk"
                    );

                    Info<< "    Writing " << currentSet.name()
                        << " (size "
                        << returnReduce(currentSet.size(), sumOp<label>())
                        << ") to "
                        << currentSet.instance()/currentSet.local()
                           /currentSet.name()
                        << " and to vtk file " << vtkName << endl << endl;

                    writeVTK(mesh, currentSet, vtkName);
                }
                else
                {
                    Info<< "    Writing " << currentSet.name()
                        << " (size "
                        << returnReduce(currentSet.size(), sumOp<label>())
                        << ") to "
                        << currentSet.instance()/currentSet.local()
                           /currentSet.name() << endl << endl;
                }

                if (writeCurrentTime)
                {
                    currentSet.instance() = mesh.time().timeName();
                }
                currentSet.write();
            }
        }
    }
    catch (Foam::IOerror& fIOErr)
    {
        ok = false;

        Pout<< fIOErr.message().c_str() << endl;

        if (sourceType.size())
        {
            Pout<< topoSetSource::usage(sourceType).c_str();
        }
    }
    catch (Foam::error& fErr)
    {
        ok = false;

        Pout<< fErr.message().c_str() << endl;

        if (sourceType.size())
        {
            Pout<< topoSetSource::usage(sourceType).c_str();
        }
    }

    return ok;
}


// Status returned from parsing the first token of the line
enum commandStatus
{
    QUIT,           // quit program
    INVALID,        // token is not a valid set manipulation command
    VALIDSETCMD,    // ,,    is a valid     ,,
    VALIDZONECMD    // ,,    is a valid     zone      ,,
};


void printMesh(const Time& runTime, const polyMesh& mesh)
{
    Info<< "Time:" << runTime.timeName()
        << "  cells:" << mesh.globalData().nTotalCells()
        << "  faces:" << mesh.globalData().nTotalFaces()
        << "  points:" << mesh.globalData().nTotalPoints()
        << "  patches:" << mesh.boundaryMesh().size()
        << "  bb:" << mesh.bounds() << nl;
}


polyMesh::readUpdateState meshReadUpdate(polyMesh& mesh)
{
    polyMesh::readUpdateState stat = mesh.readUpdate();

    switch(stat)
    {
        case polyMesh::UNCHANGED:
        {
            Info<< "    mesh not changed." << endl;
            break;
        }
        case polyMesh::POINTS_MOVED:
        {
            Info<< "    points moved; topology unchanged." << endl;
            break;
        }
        case polyMesh::TOPO_CHANGE:
        {
            Info<< "    topology changed; patches unchanged." << nl
                << "    ";
            printMesh(mesh.time(), mesh);
            break;
        }
        case polyMesh::TOPO_PATCH_CHANGE:
        {
            Info<< "    topology changed and patches changed." << nl
                << "    ";
            printMesh(mesh.time(), mesh);

            break;
        }
        default:
        {
            FatalErrorIn("meshReadUpdate(polyMesh&)")
                << "Illegal mesh update state "
                << stat  << abort(FatalError);
            break;
        }
    }
    return stat;
}


commandStatus parseType
(
    Time& runTime,
    polyMesh& mesh,
    const word& setType,
    IStringStream& is
)
{
    if (setType.empty())
    {
        Info<< "Type 'help' for usage information" << endl;

        return INVALID;
    }
    else if (setType == "help")
    {
        printHelp(Info);

        return INVALID;
    }
    else if (setType == "list")
    {
        printAllSets(mesh, Info);

        return INVALID;
    }
    else if (setType == "time")
    {
        scalar requestedTime = readScalar(is);
        instantList Times = runTime.times();

        label nearestIndex = Time::findClosestTimeIndex(Times, requestedTime);

        Info<< "Changing time from " << runTime.timeName()
            << " to " << Times[nearestIndex].name()
            << endl;

        // Set time
        runTime.setTime(Times[nearestIndex], nearestIndex);
        // Optionally re-read mesh
        meshReadUpdate(mesh);

        return INVALID;
    }
    else if (setType == "quit")
    {
        Info<< "Quitting ..." << endl;

        return QUIT;
    }
    else if
    (
        setType == "cellSet"
     || setType == "faceSet"
     || setType == "pointSet"
    )
    {
        return VALIDSETCMD;
    }
    else if
    (
        setType == "cellZoneSet"
     || setType == "faceZoneSet"
     || setType == "pointZoneSet"
    )
    {
        return VALIDZONECMD;
    }
    else
    {
        SeriousErrorIn
        (
            "commandStatus parseType(Time&, polyMesh&, const word&"
            ", IStringStream&)"
        )   << "Illegal command " << setType << endl
            << "Should be one of 'help', 'list', 'time' or a set type :"
            << " 'cellSet', 'faceSet', 'pointSet', 'faceZoneSet'"
            << endl;

        return INVALID;
    }
}


commandStatus parseAction(const word& actionName)
{
    commandStatus stat = INVALID;

    if (actionName.size())
    {
        try
        {
            (void)topoSetSource::toAction(actionName);

            stat = VALIDSETCMD;
        }
        catch (Foam::IOerror& fIOErr)
        {
            stat = INVALID;
        }
        catch (Foam::error& fErr)
        {
            stat = INVALID;
        }
    }
    return stat;
}



int main(int argc, char *argv[])
{
    timeSelector::addOptions(true, false);
#   include "addRegionOption.H"
    argList::addBoolOption("noVTK", "do not write VTK files");
    argList::addBoolOption("loop", "execute batch commands for all timesteps");
    argList::addOption
    (
        "batch",
        "file",
        "process in batch mode, using input from specified file"
    );
    argList::addBoolOption
    (
        "noSync",
        "do not synchronise selection across coupled patches"
    );

#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);

    const bool writeVTK = !args.optionFound("noVTK");
    const bool loop = args.optionFound("loop");
    const bool batch = args.optionFound("batch");
    const bool noSync = args.optionFound("noSync");

    if (loop && !batch)
    {
        FatalErrorIn(args.executable())
            << "Can only loop in batch mode."
            << exit(FatalError);
    }


#   include "createNamedPolyMesh.H"

    // Print some mesh info
    printMesh(runTime, mesh);

    // Print current sets
    printAllSets(mesh, Info);

    // Read history if interactive
#   ifdef HAS_READLINE
    if (!batch && !read_history((runTime.path()/historyFile).c_str()))
    {
        Info<< "Successfully read history from " << historyFile << endl;
    }
#   endif


    // Exit status
    int status = 0;


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        // Handle geometry/topology changes
        meshReadUpdate(mesh);


        // Main command read & execute loop
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        autoPtr<IFstream> fileStreamPtr(NULL);

        if (batch)
        {
            const fileName batchFile = args["batch"];

            Info<< "Reading commands from file " << batchFile << endl;

            // we cannot handle .gz files
            if (!isFile(batchFile, false))
            {
                FatalErrorIn(args.executable())
                    << "Cannot open file " << batchFile << exit(FatalError);
            }

            fileStreamPtr.reset(new IFstream(batchFile));
        }

        Info<< "Please type 'help', 'quit' or a set command after prompt."
            << endl;

        // Whether to quit
        bool quit = false;

        FatalError.throwExceptions();
        FatalIOError.throwExceptions();

        do
        {
            string rawLine;

            // Type: cellSet, faceSet, pointSet
            word setType;
            // Name of destination set.
            word setName;
            // Action (new, invert etc.)
            word actionName;

            commandStatus stat = INVALID;

            if (fileStreamPtr.valid())
            {
                if (!fileStreamPtr().good())
                {
                    Info<< "End of batch file" << endl;
                    // No error.
                    break;
                }

                fileStreamPtr().getLine(rawLine);

                if (rawLine.size())
                {
                    Info<< "Doing:" << rawLine << endl;
                }
            }
            else
            {
#               ifdef HAS_READLINE
                {
                    char* linePtr = readline("readline>");

                    if (linePtr)
                    {
                        rawLine = string(linePtr);

                        if (*linePtr)
                        {
                            add_history(linePtr);
                            write_history(historyFile);
                        }

                        free(linePtr);   // readline uses malloc, not new.
                    }
                    else
                    {
                        break;
                    }
                }
#               else
                {
                    if (!std::cin.good())
                    {
                        Info<< "End of cin" << endl;
                        // No error.
                        break;
                    }
                    Info<< "Command>" << flush;
                    std::getline(std::cin, rawLine);
                }
#               endif
            }

            // Strip off anything after #
            string::size_type i = rawLine.find_first_of("#");
            if (i != string::npos)
            {
                rawLine = rawLine(0, i);
            }

            if (rawLine.empty())
            {
                continue;
            }

            IStringStream is(rawLine + ' ');

            // Type: cellSet, faceSet, pointSet, faceZoneSet
            is  >> setType;

            stat = parseType(runTime, mesh, setType, is);

            if (stat == VALIDSETCMD || stat == VALIDZONECMD)
            {
                if (is >> setName)
                {
                    if (is >> actionName)
                    {
                        stat = parseAction(actionName);
                    }
                }
            }

            if (stat == QUIT)
            {
                // Make sure to quit
                quit = true;
            }
            else if (stat == VALIDSETCMD || stat == VALIDZONECMD)
            {
                bool ok = doCommand
                (
                    mesh,
                    setType,
                    setName,
                    actionName,
                    writeVTK,
                    loop,   // if in looping mode dump sets to time directory
                    noSync,
                    is
                );

                if (!ok && batch)
                {
                    // Exit with error.
                    quit = true;
                    status = 1;
                }
            }

        } while (!quit);

        if (quit)
        {
            break;
        }
    }

    Info<< "\nEnd\n" << endl;

    return status;
}


// ************************************************************************* //
