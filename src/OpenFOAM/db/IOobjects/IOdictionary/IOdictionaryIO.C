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

#include "IOdictionary.H"
#include "Pstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Parallel aware reading, using non-virtual type information (typeName instead
// of type()) because of use in constructor.
void Foam::IOdictionary::readFile(const bool masterOnly)
{
    if (Pstream::master() || !masterOnly)
    {
        if (debug)
        {
            Pout<< "IOdictionary : Reading " << objectPath()
                << " from file " << endl;
        }

        // Set flag for e.g. codeStream
        bool oldFlag = regIOobject::masterOnlyReading;
        regIOobject::masterOnlyReading = masterOnly;

        // Read file
        readStream(typeName) >> *this;
        close();

        regIOobject::masterOnlyReading = oldFlag;

        if (writeDictionaries && Pstream::master())
        {
            Sout<< nl
                << "--- IOdictionary " << name()
                << ' ' << objectPath() << ":" << nl;
            writeHeader(Sout);
            writeData(Sout);
            Sout<< "--- End of IOdictionary " << name() << nl << endl;
        }
    }

    if (masterOnly && Pstream::parRun())
    {
        // Scatter master data using communication scheme

        const List<Pstream::commsStruct>& comms =
        (
            (Pstream::nProcs() < Pstream::nProcsSimpleSum)
          ? Pstream::linearCommunication()
          : Pstream::treeCommunication()
        );

        // Master reads headerclassname from file. Make sure this gets
        // transfered as well as contents.
        Pstream::scatter
        (
            comms,
            const_cast<word&>(headerClassName()),
            Pstream::msgType(),
            Pstream::worldComm
        );
        Pstream::scatter(comms, note(), Pstream::msgType(), Pstream::worldComm);

        // Get my communication order
        const Pstream::commsStruct& myComm = comms[Pstream::myProcNo()];

        // Reveive from up
        if (myComm.above() != -1)
        {
            if (debug)
            {
                Pout<< "IOdictionary : Reading " << objectPath()
                    << " from processor " << myComm.above() << endl;
            }

            // Note: use ASCII for now - binary IO of dictionaries is
            // not currently supported or rather the primitiveEntries of
            // the dictionary think they are in binary form whereas they are
            // not. Could reset all the ITstreams to ascii?
            IPstream fromAbove
            (
                Pstream::scheduled,
                myComm.above(),
                0,
                Pstream::msgType(),
                Pstream::worldComm,
                IOstream::ASCII
            );
            IOdictionary::readData(fromAbove);
        }

        // Send to my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            if (debug)
            {
                Pout<< "IOdictionary : Sending " << objectPath()
                    << " to processor " << myComm.below()[belowI] << endl;
            }
            OPstream toBelow
            (
                Pstream::scheduled,
                myComm.below()[belowI],
                0,
                Pstream::msgType(),
                Pstream::worldComm,
                IOstream::ASCII
            );
            IOdictionary::writeData(toBelow);
        }
    }
}


// * * * * * * * * * * * * * * * Members Functions * * * * * * * * * * * * * //

bool Foam::IOdictionary::readData(Istream& is)
{
    is >> *this;

    if (writeDictionaries && Pstream::master() && !is.bad())
    {
        Sout<< nl
            << "--- IOdictionary " << name()
            << ' ' << objectPath() << ":" << nl;
        writeHeader(Sout);
        writeData(Sout);
        Sout<< "--- End of IOdictionary " << name() << nl << endl;
    }

    return !is.bad();
}


bool Foam::IOdictionary::writeData(Ostream& os) const
{
    dictionary::write(os, false);
    return os.good();
}


// ************************************************************************* //
