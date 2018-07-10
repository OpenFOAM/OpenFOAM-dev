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

Description
    Create intermediate mesh from PROSTAR files

\*---------------------------------------------------------------------------*/

#include "starMesh.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::starMesh::readCouples()
{
    fileName couplesFileName(casePrefix_ + ".cpl");

    label nCouples = 0;

    // Count number of couples
    {
        IFstream couplesFile(couplesFileName);

        if (couplesFile.good())
        {
            Info<< "\nReading couples" << endl;

            label matchLabel, nEntries, typeFlag;
            label starMasterCell, rotXMasterFace;
            label starSlaveCell, rotXSlaveFace;

            // count the number of entries to read
            while (!(couplesFile >> matchLabel).eof())
            {
                // read number of entries and match type flag.
                couplesFile >> nEntries;

                couplesFile >> typeFlag;

                // read master cell and face
                couplesFile >> starMasterCell >> rotXMasterFace;

                // add number of couples from current match
                label nSlavesToRead = nEntries - 1;

                nCouples += nSlavesToRead;

                for (int i = 0; i < nSlavesToRead; i++)
                {
                    couplesFile >> starSlaveCell >> rotXSlaveFace;
                }
            }

            Info<< "Number of couples = " << nCouples << endl << endl;
        }
        else
        {
            Info<< endl << "No couple matches defined." << endl;
        }
    }

    // Read couples
    if (nCouples > 0)
    {
        // read couples
        couples_.setSize(nCouples);
        label couplei = 0;

        // A mesh with couples cannot be a shape mesh
        isShapeMesh_ = false;

        IFstream couplesFile(couplesFileName);

        label matchLabel, nEntries, typeFlag;
        label starMasterCell, masterCell, rotXMasterFace, rotZeroMasterFace;
        label starSlaveCell, slaveCell, rotXSlaveFace, rotZeroSlaveFace;

        while (!(couplesFile >> matchLabel).eof())
        {
            // read number of entries and match type flag.
            // Note. At the moment, only integral matches are supported
            couplesFile >> nEntries;

            couplesFile >> typeFlag;

            // read master cell and face
            couplesFile >> starMasterCell >> rotXMasterFace;

            // translate the cell labels
            masterCell = starCellLabelLookup_[starMasterCell];

            // translate the master face into rotation zero if applicable
            if (starCellPermutation_[masterCell] > -1)
            {
                const label curMasterPermutation =
                    starCellPermutation_[masterCell];

                rotZeroMasterFace =
                    sammFacePermutationTable
                        [curMasterPermutation]
                        [rotXMasterFace];
            }
            else
            {
                rotZeroMasterFace = rotXMasterFace;
            }

            // get master face index
            label masterFaceID =
                shapeFaceLookup
                    [cellShapes_[masterCell].model().index()]
                    [rotZeroMasterFace];

            // number of slave faces
            label nSlavesToRead = nEntries - 1;

            for (int i = 0; i < nSlavesToRead; i++)
            {
                couplesFile >> starSlaveCell >> rotXSlaveFace;

                // translate the cell labels
                slaveCell = starCellLabelLookup_[starSlaveCell];

                // translate the slave face into rotation zero if applicable
                if (starCellPermutation_[slaveCell] > -1)
                {
                    const label curSlavePermutation =
                        starCellPermutation_[slaveCell];

                    rotZeroSlaveFace =
                        sammFacePermutationTable
                            [curSlavePermutation]
                            [rotXSlaveFace];
                }
                else
                {
                    rotZeroSlaveFace = rotXSlaveFace;
                }

                label slaveFaceID =
                    shapeFaceLookup
                        [cellShapes_[slaveCell].model().index()]
                        [rotZeroSlaveFace];

                // Set the couple
                couples_.set
                (
                    couplei++,
                    new coupledFacePair
                    (
                        matchLabel,
                        masterCell, masterFaceID,
                        slaveCell, slaveFaceID,
                        typeFlag
                    )
                );
            }
        }

        Info<< "finished reading couples" << endl;
    }
}


// ************************************************************************* //
