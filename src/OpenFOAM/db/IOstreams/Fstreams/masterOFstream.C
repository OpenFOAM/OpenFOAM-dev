/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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

#include "masterOFstream.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "PstreamBuffers.H"
#include "masterUncollatedFileOperation.H"
#include "boolList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::masterOFstream::checkWrite
(
    const fileName& fName,
    const string& str
)
{
    mkDir(fName.path());

    OFstream os
    (
        fName,
        IOstream::BINARY,
        version(),
        compression_,
        append_
    );

    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Could not open file " << fName
            << exit(FatalIOError);
    }

    os.writeQuoted(str, false);
    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Failed writing to " << fName
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::masterOFstream::masterOFstream
(
    const fileName& filePath,
    streamFormat format,
    versionNumber version,
    compressionType compression,
    const bool append,
    const bool write
)
:
    OStringStream(format, version),
    filePath_(filePath),
    compression_(compression),
    append_(append),
    write_(write)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::masterOFstream::~masterOFstream()
{
    if (Pstream::parRun())
    {
        List<fileName> filePaths(Pstream::nProcs());
        filePaths[Pstream::myProcNo()] = filePath_;
        Pstream::gatherList(filePaths);

        bool uniform =
            fileOperations::masterUncollatedFileOperation::uniformFile
            (
                filePaths
            );

        Pstream::scatter(uniform);

        if (uniform)
        {
            if (Pstream::master() && write_)
            {
                checkWrite(filePath_, str());
            }
            return;
        }
        boolList write(Pstream::nProcs());
        write[Pstream::myProcNo()] = write_;
        Pstream::gatherList(write);


        // Different files
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send my buffer to master
        if (!Pstream::master())
        {
            UOPstream os(Pstream::masterNo(), pBufs);
            string s(this->str());
            os.write(&s[0], s.size());
        }

        labelList recvSizes;
        pBufs.finishedSends(recvSizes);

        if (Pstream::master())
        {
            // Write my own data
            {
                if (write[Pstream::myProcNo()])
                {
                    checkWrite(filePaths[Pstream::myProcNo()], str());
                }
            }

            for (label proci = 1; proci < Pstream::nProcs(); proci++)
            {
                UIPstream is(proci, pBufs);
                List<char> buf(recvSizes[proci]);

                is.read(buf.begin(), buf.size());

                if (write[proci])
                {
                    checkWrite
                    (
                        filePaths[proci],
                        string(buf.begin(), buf.size())
                    );
                }
            }
        }
    }
    else
    {
        checkWrite(filePath_, str());
    }
}


// ************************************************************************* //
