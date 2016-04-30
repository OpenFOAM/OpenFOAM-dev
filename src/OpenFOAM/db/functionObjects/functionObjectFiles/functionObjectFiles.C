/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "functionObjectFiles.H"
#include "Time.H"
#include "IFstream.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjectFiles::createFiles()
{
    if (Pstream::master())
    {
        const word startTimeName =
            obr_.time().timeName(obr_.time().startTime().value());

        forAll(names_, i)
        {
            if (!filePtrs_.set(i))
            {
                fileName outputDir(baseFileDir()/prefix_/startTimeName);
                mkDir(outputDir);

                word fName(names_[i]);

                // Check if file already exists
                IFstream is(outputDir/(fName + ".dat"));
                if (is.good())
                {
                    fName = fName + "_" + obr_.time().timeName();
                }

                filePtrs_.set(i, new OFstream(outputDir/(fName + ".dat")));

                initStream(filePtrs_[i]);

                writeFileHeader(i);

            }
        }
    }
}


void Foam::functionObjectFiles::write()
{
    createFiles();
}


void Foam::functionObjectFiles::resetNames(const wordList& names)
{
    names_.clear();
    names_.append(names);

    if (Pstream::master())
    {
        filePtrs_.clear();
        filePtrs_.setSize(names_.size());

        createFiles();
    }
}


void Foam::functionObjectFiles::resetName(const word& name)
{
    names_.clear();
    names_.append(name);

    if (Pstream::master())
    {
        filePtrs_.clear();
        filePtrs_.setSize(1);

        createFiles();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjectFiles::functionObjectFiles
(
    const objectRegistry& obr,
    const word& prefix
)
:
    functionObjectFile(obr, prefix),
    names_(),
    filePtrs_()
{}


Foam::functionObjectFiles::functionObjectFiles
(
    const objectRegistry& obr,
    const word& prefix,
    const word& name
)
:
    functionObjectFile(obr, prefix),
    names_(),
    filePtrs_()
{
    names_.clear();
    names_.append(name);
    if (Pstream::master())
    {
        filePtrs_.clear();
        filePtrs_.setSize(1);

        // Cannot create files - need to access virtual function
    }
}


Foam::functionObjectFiles::functionObjectFiles
(
    const objectRegistry& obr,
    const word& prefix,
    const wordList& names
)
:
    functionObjectFile(obr, prefix),
    names_(names),
    filePtrs_()
{
    names_.clear();
    names_.append(names);
    if (Pstream::master())
    {
        filePtrs_.clear();
        filePtrs_.setSize(names_.size());

        // Cannot create files - need to access virtual function
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjectFiles::~functionObjectFiles()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::functionObjectFiles::names() const
{
    return names_;
}


Foam::OFstream& Foam::functionObjectFiles::file()
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


Foam::PtrList<Foam::OFstream>& Foam::functionObjectFiles::files()
{
    if (!Pstream::master())
    {
        FatalErrorInFunction
            << "Request for files() can only be done by the master process"
            << abort(FatalError);
    }

    return filePtrs_;
}


Foam::OFstream& Foam::functionObjectFiles::file(const label i)
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


// ************************************************************************* //
