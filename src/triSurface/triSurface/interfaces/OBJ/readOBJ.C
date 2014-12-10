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

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "IFstream.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::triSurface::readOBJ(const fileName& OBJfileName)
{
    IFstream OBJfile(OBJfileName);

    if (!OBJfile.good())
    {
        FatalErrorIn("triSurface::readOBJ(const fileName&)")
            << "Cannot read file " << OBJfileName
            << exit(FatalError);
    }

    DynamicList<point> points;
    DynamicList<labelledTri> faces;
    HashTable<label> groupToPatch;

    label groupID = 0;
    label maxGroupID = 0;

    while (OBJfile.good())
    {
        string line = getLineNoComment(OBJfile);

        label sz = line.size();

        if (sz && line[sz-1] == '\\')
        {
            line.substr(0, sz-1);
            line += getLineNoComment(OBJfile);
        }

        // Read first word
        IStringStream lineStream(line);
        word cmd;
        lineStream >> cmd;

        if (cmd == "v")
        {
            scalar x, y, z;

            lineStream >> x >> y >> z;

            points.append(point(x, y, z));
        }
        else if (cmd == "g")
        {
            word group;

            lineStream >> group;

            HashTable<label>::const_iterator findGroup =
                groupToPatch.find(group);

            if (findGroup != groupToPatch.end())
            {
                groupID = findGroup();
            }
            else
            {
                groupID = maxGroupID;

                groupToPatch.insert(group, groupID);

                maxGroupID++;
            }
        }
        else if (cmd == "f")
        {
            DynamicList<label> verts;

            // Assume 'f' is followed by space.
            string::size_type endNum = 1;

            while (true)
            {
                string::size_type startNum =
                    line.find_first_not_of(" \r", endNum);

                if (startNum == string::npos)
                {
                    break;
                }

                endNum = line.find(' ', startNum);

                string vertexSpec;
                if (endNum != string::npos)
                {
                    vertexSpec = line.substr(startNum, endNum-startNum);
                }
                else
                {
                    vertexSpec = line.substr(startNum, line.size() - startNum);
                }

                string::size_type slashPos = vertexSpec.find('/');

                label vertI = 0;
                if (slashPos != string::npos)
                {
                    IStringStream intStream(vertexSpec.substr(0, slashPos));

                    intStream >> vertI;
                }
                else
                {
                    IStringStream intStream(vertexSpec);

                    intStream >> vertI;
                }
                verts.append(vertI - 1);
            }

            verts.shrink();

            // Do simple face triangulation around f[0].
            // Cannot use face::triangulation since no complete points yet.
            for (label fp = 1; fp < verts.size() - 1; fp++)
            {
                label fp1 = verts.fcIndex(fp);

                labelledTri tri(verts[0], verts[fp], verts[fp1], groupID);

                faces.append(tri);
            }
        }
    }

    points.shrink();
    faces.shrink();

    // Convert groupToPatch to patchList.
    geometricSurfacePatchList patches(maxGroupID);

    if (maxGroupID == 0)
    {
        // Generate default patch
        patches.setSize(1);
        patches[0] = geometricSurfacePatch("empty", "patch0", 0);
    }
    else
    {
        forAllConstIter(HashTable<label>, groupToPatch, iter)
        {
            patches[iter()] = geometricSurfacePatch
            (
                "empty",
                iter.key(),
                iter()
            );
        }
    }


    // Transfer DynamicLists to straight ones.
    pointField allPoints(points.xfer());

    // Create triSurface
    *this = triSurface(faces, patches, allPoints, true);

    return true;
}


// ************************************************************************* //
