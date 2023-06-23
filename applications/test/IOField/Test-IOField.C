/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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
    Test-IOField

Description
    Test the processor-local reading of IOField (used in the lagrangian libs)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "IOField.H"
#include "polyMesh.H"
#include "Time.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void write(const typeIOobject<IOField<label>>& io, const label sz)
{
    IOField<label> fld(io, sz);
    forAll(fld, i)
    {
        fld[i] = i+1000;
    }
    Pout<< "writing:" << fld << endl;
    fld.write(sz > 0);
}


void read(const typeIOobject<IOField<label>>& io, const label sz)
{
    bool valid = (sz > 0);
    Pout<< "    valid:" << valid << endl;
    IOField<label> fld(io, valid);
    Pout<< "    wanted:" << sz << " actually read:" << fld.size() << endl;

    if (fld.size() != sz)
    {
        FatalErrorInFunction<< "io:" << fld.objectPath() << exit(FatalError);
    }
}


void writeAndRead
(
    const typeIOobject<IOField<label>>& io,
    const label sz,
    const word& writeType,
    const IOobject::readOption readOpt,
    const word& readType
)
{
    Pout<< "** Writing:" << writeType
        << " Reading:" << readType << endl;

    autoPtr<fileOperation> writeHandler
    (
        fileOperation::New(writeType, writeInfoHeader)
    );
    fileHandler(writeHandler);

    // Delete
    Pout<< "Deleting:" << io.filePath() << endl;
    fileHandler().rm(io.objectPath());

    // Write
    Pout<< "Writing:" << io.objectPath() << endl;
    write(io, sz);

    autoPtr<fileOperation> readHandler
    (
        fileOperation::New(readType, writeInfoHeader)
    );
    fileHandler(readHandler);

    // Read
    typeIOobject<IOField<label>> readIO(io);
    readIO.readOpt() = readOpt;
    Pout<< "Reading:" << readIO.filePath() << endl;
    read(readIO, sz);

    Pout<< "** Done writing:" << writeType
        << " Reading:" << readType << endl << endl << endl;
}


void readIfPresent
(
    typeIOobject<IOField<label>>& io,
    const label sz,
    const word& readType
)
{
    autoPtr<fileOperation> readHandler
    (
        fileOperation::New(readType, writeInfoHeader)
    );
    fileHandler(readHandler);

    // Read
    Pout<< "Reading:" << io.filePath() << endl;
    io.readOpt() = IOobject::READ_IF_PRESENT;
    read(io, sz);
}


// Main program:

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"
    timeSelector::select0(runTime, args);
    #include "createPolyMesh.H"

    label sz = 0;
    if (Pstream::myProcNo() % 2)
    {
        sz = 1;
    }

    typeIOobject<IOField<label>> io
    (
        "bla",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    wordList handlers
    (
        Foam::fileOperation::wordConstructorTablePtr_->sortedToc()
    );

    Info<< "Found handlers:" << handlers << endl;

    forAll(handlers, writei)
    {
        const word& writeHandler = handlers[writei];

        forAll(handlers, readi)
        {
            const word& readHandler = handlers[readi];

            writeAndRead
            (
                io,
                sz,
                writeHandler,
                IOobject::READ_IF_PRESENT,
                readHandler
            );
        }
    }


    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
