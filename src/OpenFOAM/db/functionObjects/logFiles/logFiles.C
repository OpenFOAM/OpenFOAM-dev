/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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

#include "logFiles.H"
#include "Time.H"
#include "OSspecific.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::logFiles::createFiles()
{
    if (Pstream::master())
    {
        const word timeName = fileObr_.time().timeName();

        forAll(names_, i)
        {
            if (!filePtrs_.set(i))
            {
                const fileName outputDir(baseFileDir()/prefix_/timeName);
                mkDir(outputDir);
                filePtrs_.set(i, new OFstream(outputDir/(names_[i] + ".dat")));
                initStream(filePtrs_[i]);
                writeFileHeader(i);
            }
        }
    }
}


void Foam::functionObjects::logFiles::resetNames(const wordList& names)
{
    names_.clear();
    names_.append(names);

    if (Pstream::master())
    {
        filePtrs_.clear();
        filePtrs_.setSize(names_.size());
    }
}


void Foam::functionObjects::logFiles::resetName(const word& name)
{
    names_.clear();
    names_.append(name);

    if (Pstream::master())
    {
        filePtrs_.clear();
        filePtrs_.setSize(1);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::logFiles::logFiles
(
    const objectRegistry& obr,
    const word& prefix
)
:
    writeFile(obr, prefix),
    names_(),
    filePtrs_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::logFiles::~logFiles()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::functionObjects::logFiles::names() const
{
    return names_;
}


Foam::OFstream& Foam::functionObjects::logFiles::file()
{
    if (!Pstream::master())
    {
        FatalErrorInFunction
            << "Request for file() can only be done by the master process"
            << abort(FatalError);
    }

    if (filePtrs_.size() != 1)
    {
        WarningInFunction
            << "Requested single file, but multiple files are present"
            << endl;
    }

    if (!filePtrs_.set(0))
    {
        FatalErrorInFunction
            << "File pointer at index " << 0 << " not allocated"
            << abort(FatalError);
    }

    return filePtrs_[0];
}


Foam::PtrList<Foam::OFstream>& Foam::functionObjects::logFiles::files()
{
    if (!Pstream::master())
    {
        FatalErrorInFunction
            << "Request for files() can only be done by the master process"
            << abort(FatalError);
    }

    return filePtrs_;
}


Foam::OFstream& Foam::functionObjects::logFiles::file(const label i)
{
    if (!Pstream::master())
    {
        FatalErrorInFunction
            << "Request for file(i) can only be done by the master process"
            << abort(FatalError);
    }

    if (!filePtrs_.set(i))
    {
        FatalErrorInFunction
            << "File pointer at index " << i << " not allocated"
            << abort(FatalError);
    }

    return filePtrs_[i];
}


bool Foam::functionObjects::logFiles::write()
{
    createFiles();

    return true;
}


// ************************************************************************* //
