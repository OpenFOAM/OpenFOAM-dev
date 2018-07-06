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

#include "writeFaceSet.H"
#include "OFstream.H"
#include "writeFuns.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::writeFaceSet
(
    const bool binary,
    const vtkMesh& vMesh,
    const faceSet& set,
    const fileName& fileName
)
{
    const faceList& faces = vMesh.mesh().faces();

    std::ofstream ostr(fileName.c_str());

    writeFuns::writeHeader
    (
        ostr,
        binary,
        set.name()
    );

    ostr<< "DATASET POLYDATA" << std::endl;

    //------------------------------------------------------------------
    //
    // Write topology
    //
    //------------------------------------------------------------------


    // Construct primitivePatch of faces in faceSet.

    faceList setFaces(set.size());
    labelList setFaceLabels(set.size());
    label setFacei = 0;

    forAllConstIter(faceSet, set, iter)
    {
        setFaceLabels[setFacei] = iter.key();
        setFaces[setFacei] = faces[iter.key()];
        setFacei++;
    }
    primitiveFacePatch fp(setFaces, vMesh.mesh().points());


    // Write points and faces as polygons

    ostr<< "POINTS " << fp.nPoints() << " float" << std::endl;

    DynamicList<floatScalar> ptField(3*fp.nPoints());

    writeFuns::insert(fp.localPoints(), ptField);

    writeFuns::write(ostr, binary, ptField);


    label nFaceVerts = 0;

    forAll(fp.localFaces(), facei)
    {
        nFaceVerts += fp.localFaces()[facei].size() + 1;
    }
    ostr<< "POLYGONS " << fp.size() << ' ' << nFaceVerts << std::endl;


    DynamicList<label> vertLabels(nFaceVerts);

    forAll(fp.localFaces(), facei)
    {
        const face& f = fp.localFaces()[facei];

        vertLabels.append(f.size());

        writeFuns::insert(f, vertLabels);
    }
    writeFuns::write(ostr, binary, vertLabels);


    //-----------------------------------------------------------------
    //
    // Write data
    //
    //-----------------------------------------------------------------

    // Write faceID

    ostr
        << "CELL_DATA " << fp.size() << std::endl
        << "FIELD attributes 1" << std::endl;

    // Cell ids first
    ostr<< "faceID 1 " << fp.size() << " int" << std::endl;

    writeFuns::write(ostr, binary, setFaceLabels);
}


// ************************************************************************* //
