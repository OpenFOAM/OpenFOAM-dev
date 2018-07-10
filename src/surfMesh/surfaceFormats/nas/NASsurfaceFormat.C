/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "NASsurfaceFormat.H"
#include "IFstream.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::NASsurfaceFormat<Face>::NASsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::NASsurfaceFormat<Face>::read
(
    const fileName& filename
)
{
    const bool mustTriangulate = this->isTri();
    this->clear();

    IFstream is(filename);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << filename
            << exit(FatalError);
    }

    // Nastran index of points
    DynamicList<label>  pointId;
    DynamicList<point>  dynPoints;
    DynamicList<Face>   dynFaces;
    DynamicList<label>  dynZones;
    DynamicList<label>  dynSizes;
    Map<label>          lookup;

    // assume the types are not intermixed
    // leave faces that didn't have a group in 0
    bool sorted = true;
    label zoneI = 0;

    // Name for face group
    Map<word> nameLookup;

    // Ansa tags. Denoted by $ANSA_NAME.
    // These will appear just before the first use of a type.
    // We read them and store the PSHELL types which are used to name
    // the zones.
    label ansaId = -1;
    word  ansaType, ansaName;

    // A single warning per unrecognized command
    HashSet<word> unhandledCmd;

    while (is.good())
    {
        string line;
        is.getLine(line);

        // Ansa extension
        if (line.substr(0, 10) == "$ANSA_NAME")
        {
            string::size_type sem0 = line.find (';', 0);
            string::size_type sem1 = line.find (';', sem0+1);
            string::size_type sem2 = line.find (';', sem1+1);

            if
            (
                sem0 != string::npos
             && sem1 != string::npos
             && sem2 != string::npos
            )
            {
                ansaId = readLabel
                (
                    IStringStream(line.substr(sem0+1, sem1-sem0-1))()
                );
                ansaType = line.substr(sem1+1, sem2-sem1-1);

                string rawName;
                is.getLine(rawName);
                if (rawName[rawName.size()-1] == '\r')
                {
                    rawName = rawName.substr(1, rawName.size()-2);
                }
                else
                {
                    rawName = rawName.substr(1, rawName.size()-1);
                }

                string::stripInvalid<word>(rawName);
                ansaName = rawName;

                // Info<< "ANSA tag for NastranID:" << ansaId
                //     << " of type " << ansaType
                //     << " name " << ansaName << endl;
            }
        }


        // Hypermesh extension
        // $HMNAME COMP                   1"partName"
        if
        (
            line.substr(0, 12) == "$HMNAME COMP"
         && line.find ('"') != string::npos
        )
        {
            label groupId = readLabel
            (
                IStringStream(line.substr(16, 16))()
            );

            IStringStream lineStream(line.substr(32));

            string rawName;
            lineStream >> rawName;
            string::stripInvalid<word>(rawName);

            word groupName(rawName);
            nameLookup.insert(groupId, groupName);

            // Info<< "group " << groupId << " => " << groupName << endl;
        }


        // Skip empty or comment
        if (line.empty() || line[0] == '$')
        {
            continue;
        }

        // Check if character 72 is continuation
        if (line.size() > 72 && line[72] == '+')
        {
            line = line.substr(0, 72);

            while (true)
            {
                string buf;
                is.getLine(buf);

                if (buf.size() > 72 && buf[72] == '+')
                {
                    line += buf.substr(8, 64);
                }
                else
                {
                    line += buf.substr(8, buf.size()-8);
                    break;
                }
            }
        }


        // Read first word
        IStringStream lineStream(line);
        word cmd;
        lineStream >> cmd;

        if (cmd == "CTRIA3")
        {
            triFace fTri;

            label groupId = readLabel(IStringStream(line.substr(16,8))());
            fTri[0] = readLabel(IStringStream(line.substr(24,8))());
            fTri[1] = readLabel(IStringStream(line.substr(32,8))());
            fTri[2] = readLabel(IStringStream(line.substr(40,8))());

            // Convert groupID into zoneId
            Map<label>::const_iterator fnd = lookup.find(groupId);
            if (fnd != lookup.end())
            {
                if (zoneI != fnd())
                {
                    // pshell types are intermixed
                    sorted = false;
                }
                zoneI = fnd();
            }
            else
            {
                zoneI = dynSizes.size();
                lookup.insert(groupId, zoneI);
                dynSizes.append(0);
                // Info<< "zone" << zoneI << " => group " << groupId <<endl;
            }

            dynFaces.append(fTri);
            dynZones.append(zoneI);
            dynSizes[zoneI]++;
        }
        else if (cmd == "CQUAD4")
        {
            face fQuad(4);
            labelUList& f = static_cast<labelUList&>(fQuad);

            label groupId = readLabel(IStringStream(line.substr(16,8))());
            fQuad[0] = readLabel(IStringStream(line.substr(24,8))());
            fQuad[1] = readLabel(IStringStream(line.substr(32,8))());
            fQuad[2] = readLabel(IStringStream(line.substr(40,8))());
            fQuad[3] = readLabel(IStringStream(line.substr(48,8))());

            // Convert groupID into zoneId
            Map<label>::const_iterator fnd = lookup.find(groupId);
            if (fnd != lookup.end())
            {
                if (zoneI != fnd())
                {
                    // pshell types are intermixed
                    sorted = false;
                }
                zoneI = fnd();
            }
            else
            {
                zoneI = dynSizes.size();
                lookup.insert(groupId, zoneI);
                dynSizes.append(0);
                // Info<< "zone" << zoneI << " => group " << groupId <<endl;
            }


            if (mustTriangulate)
            {
                dynFaces.append(triFace(f[0], f[1], f[2]));
                dynFaces.append(triFace(f[0], f[2], f[3]));
                dynZones.append(zoneI);
                dynZones.append(zoneI);
                dynSizes[zoneI] += 2;
            }
            else
            {
                dynFaces.append(Face(f));
                dynZones.append(zoneI);
                dynSizes[zoneI]++;
            }
        }
        else if (cmd == "GRID")
        {
            label index = readLabel(IStringStream(line.substr(8,8))());
            scalar x = parseNASCoord(line.substr(24, 8));
            scalar y = parseNASCoord(line.substr(32, 8));
            scalar z = parseNASCoord(line.substr(40, 8));

            pointId.append(index);
            dynPoints.append(point(x, y, z));
        }
        else if (cmd == "GRID*")
        {
            // Long format is on two lines with '*' continuation symbol
            // on start of second line.
            // Typical line (spaces compacted)
            // GRID*      126   0 -5.55999875E+02 -5.68730474E+02
            // *         2.14897901E+02

            label index = readLabel(IStringStream(line.substr(8,16))());
            scalar x = parseNASCoord(line.substr(40, 16));
            scalar y = parseNASCoord(line.substr(56, 16));

            is.getLine(line);
            if (line[0] != '*')
            {
                FatalErrorInFunction
                    << "Expected continuation symbol '*' when reading GRID*"
                    << " (double precision coordinate) format" << nl
                    << "Read:" << line << nl
                    << "File:" << is.name() << " line:" << is.lineNumber()
                    << exit(FatalError);
            }
            scalar z = parseNASCoord(line.substr(8, 16));

            pointId.append(index);
            dynPoints.append(point(x, y, z));
        }
        else if (cmd == "PSHELL")
        {
            // pshell type for zone names with the Ansa extension
            label groupId = readLabel(IStringStream(line.substr(8,8))());

            if (groupId == ansaId && ansaType == "PSHELL")
            {
                nameLookup.insert(ansaId, ansaName);
                // Info<< "group " << groupId << " => " << ansaName << endl;
            }
        }
        else if (unhandledCmd.insert(cmd))
        {
            Info<< "Unhandled Nastran command " << line << nl
                << "File:" << is.name() << " line:" << is.lineNumber()
                << endl;
        }
    }

    //    Info<< "Read faces:" << dynFaces.size()
    //        << " points:" << dynPoints.size()
    //        << endl;

    // transfer to normal lists
    this->storedPoints().transfer(dynPoints);

    pointId.shrink();
    dynFaces.shrink();

    // Build inverse mapping (NASTRAN pointId -> index)
    Map<label> mapPointId(2*pointId.size());
    forAll(pointId, i)
    {
        mapPointId.insert(pointId[i], i);
    }

    // Relabel faces
    // ~~~~~~~~~~~~~
    forAll(dynFaces, i)
    {
        Face& f = dynFaces[i];
        forAll(f, fp)
        {
            f[fp] = mapPointId[f[fp]];
        }
    }
    pointId.clearStorage();
    mapPointId.clear();


    // create default zone names, or from ANSA/Hypermesh information
    List<word> names(dynSizes.size());
    forAllConstIter(Map<label>, lookup, iter)
    {
        const label zoneI  = iter();
        const label groupI = iter.key();

        Map<word>::const_iterator fnd = nameLookup.find(groupI);
        if (fnd != nameLookup.end())
        {
            names[zoneI] = fnd();
        }
        else
        {
            names[zoneI] = word("zone") + ::Foam::name(zoneI);
        }
    }

    this->sortFacesAndStore(dynFaces.xfer(), dynZones.xfer(), sorted);

    // add zones, culling empty ones
    this->addZones(dynSizes, names, true);

    return true;
}


// ************************************************************************* //
