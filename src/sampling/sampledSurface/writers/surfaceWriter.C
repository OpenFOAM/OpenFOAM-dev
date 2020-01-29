/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "surfaceWriter.H"
#include "MeshedSurfaceProxy.H"
#include "proxySurfaceWriter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceWriter, 0);
    defineRunTimeSelectionTable(surfaceWriter, word);
    defineRunTimeSelectionTable(surfaceWriter, wordDict);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceWriter>
Foam::surfaceWriter::New
(
    const word& writeType,
    const IOstream::streamFormat format
)
{
    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(writeType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        if (MeshedSurfaceProxy<face>::canWriteType(writeType))
        {
            // generally unknown, but can be written via MeshedSurfaceProxy
            // use 'proxy' handler instead
            return autoPtr<surfaceWriter>(new proxySurfaceWriter(writeType));
        }

        if (cstrIter == wordConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown write type \"" << writeType << "\"\n\n"
                << "Valid write types : "
                << wordConstructorTablePtr_->sortedToc() << nl
                << "Valid proxy types : "
                << MeshedSurfaceProxy<face>::writeTypes() << endl
                << exit(FatalError);
        }
    }

    return autoPtr<surfaceWriter>(cstrIter()(format));
}


Foam::autoPtr<Foam::surfaceWriter>
Foam::surfaceWriter::New
(
    const word& writeType,
    const dictionary& optDict
)
{
    // find constructors with dictionary options
    wordDictConstructorTable::iterator cstrIter =
        wordDictConstructorTablePtr_->find(writeType);

    if (cstrIter == wordDictConstructorTablePtr_->end())
    {
        IOstream::streamFormat writeFormat = IOstream::ASCII;

        if (optDict.found("writeFormat"))
        {
            writeFormat = IOstream::formatEnum
            (
                optDict.lookup("writeFormat")
            );
        }

        // Revert to versions without options
        return Foam::surfaceWriter::New(writeType, writeFormat);
    }

    return autoPtr<surfaceWriter>(cstrIter()(optDict));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriter::surfaceWriter(const IOstream::streamFormat writeFormat)
:
    writeFormat_(writeFormat)
{}


Foam::surfaceWriter::surfaceWriter(const dictionary& optDict)
:
    writeFormat_(IOstream::ASCII)
{
    if (optDict.found("writeFormat"))
    {
        writeFormat_ = IOstream::formatEnum
        (
            optDict.lookup("writeFormat")
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceWriter::~surfaceWriter()
{}


// ************************************************************************* //
