/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description
    Create intermediate mesh from SAMM files

\*---------------------------------------------------------------------------*/

#include "sammMesh.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void sammMesh::readCouples()
{
    fileName couplesFileName(casePrefix_ + ".cpl");

    IFstream couplesFile(couplesFileName);

    if (couplesFile.good())
    {
        Info<< "\nReading couples" << endl;

        // A mesh with couples cannot be a shape mesh
        isShapeMesh_ = false;

        label matchLabel, nEntries, typeFlag;
        label masterCell, masterFace;
        label slaveCell, slaveFace;

        while (!(couplesFile >> matchLabel).eof())
        {
            // read number of entries and match type.
            // Note. At the moment, only integral matches are supported
            couplesFile >> nEntries;

            couplesFile >> typeFlag;

            if (typeFlag > 1)
            {
                Info
                    << "void sammMesh::readCouples() : "
                    << "couple " << matchLabel << " is not an integral match. "
                    << "Currently not supported" << endl;
            }

            // read master cell and face
            couplesFile >> masterCell >> masterFace;

            // get reference to master cell faces
            faceList& masterFaces = cellFaces_[masterCell - 1];

//             Info<< "Master cell: " << masterCell - 1 << " index: "
//                 << cellShapes_[masterCell - 1].model().index()
//                 << " face: " <<
//                 masterFaces
//                 [
//                     shapeFaceLookup
//                        [cellShapes_[masterCell - 1].model().index()]
//                        [masterFace]
//                 ]
//                 << endl;

            // reset master face to zero size. It cannot be removed at this
            // stage because thisw would mess up the numbering in case of
            // more than one couple an a single master cell
            masterFaces
                [
                    shapeFaceLookup
                       [cellShapes_[masterCell - 1].model().index()]
                       [masterFace]
                ].setSize(0);

            // number of slave faces
            label nSlavesToRead = nEntries - 1;

            // get index for slave face add
            label slaveToAdd = masterFaces.size();

            // reset size of master faces to accept new (couple) faces
            masterFaces.setSize(masterFaces.size() + nSlavesToRead);

            for (int i = 0; i < nSlavesToRead; i++)
            {
                couplesFile >> slaveCell >> slaveFace;

                masterFaces[slaveToAdd] =
                    cellFaces_
                        [
                            slaveCell - 1
                        ]
                        [
                            shapeFaceLookup
                                [cellShapes_[slaveCell - 1].model().index()]
                                [slaveFace]
                        ].reverseFace();

//                 Info<< " slave cell: " << slaveCell - 1 << " index: "
//                     << cellShapes_[slaveCell - 1].model().index()
//                     << " face: " << masterFaces[slaveToAdd] << endl;

                slaveToAdd++;

            }
//             Info<< endl;

        }

        // Once all couples are read, remove zero size faces from all cells
        forAll(cellFaces_, cellI)
        {
            faceList& curFaces = cellFaces_[cellI];

            label zeroSizeFound = 0;

            forAll(curFaces, faceI)
            {
                if (curFaces[faceI].empty())
                {
                    zeroSizeFound++;
                }
            }

            if (zeroSizeFound > 0)
            {
                // compress the list. A copy needs to made first
                faceList oldFaces = curFaces;

                curFaces.setSize(curFaces.size() - zeroSizeFound);

                label nFaces = 0;

                forAll(oldFaces, faceI)
                {
                    if (oldFaces[faceI].size())
                    {
                        curFaces[nFaces] = oldFaces[faceI];

                        nFaces++;
                    }
                }
            }
        }
    }
    else
    {
        Info
            << "void sammMesh::readCouples() : "
            << "Cannot read file "
            << couplesFileName
            << ". No matches defined."
            << endl;
        }
}


// ************************************************************************* //
