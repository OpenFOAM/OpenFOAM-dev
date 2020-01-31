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

#include "dlLibraryTable.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dlLibraryTable, 0);
}

Foam::dlLibraryTable Foam::libs;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dlLibraryTable::dlLibraryTable()
{}


Foam::dlLibraryTable::dlLibraryTable
(
    const dictionary& dict,
    const word& libsEntry
)
{
    open(dict, libsEntry);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dlLibraryTable::~dlLibraryTable()
{
    forAllReverse(libPtrs_, i)
    {
        if (libPtrs_[i])
        {
            if (debug)
            {
                InfoInFunction
                    << "Closing " << libNames_[i]
                    << " with handle " << uintptr_t(libPtrs_[i]) << endl;
            }
            if (!dlClose(libPtrs_[i]))
            {
                WarningInFunction<< "Failed closing " << libNames_[i]
                    << " with handle " << uintptr_t(libPtrs_[i]) << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dlLibraryTable::open
(
    const fileName& libName,
    const bool verbose
)
{
    if (libName.size())
    {
        void* libPtr = dlOpen
        (
            fileName(libName).expand(),
            verbose
        );

        if (debug)
        {
            InfoInFunction
                << "Opened " << libName
                << " resulting in handle " << uintptr_t(libPtr) << endl;
        }

        if (!libPtr)
        {
            if (verbose)
            {
                WarningInFunction
                    << "could not load " << libName
                    << endl;
            }

            return false;
        }
        else
        {
            libPtrs_.append(libPtr);
            libNames_.append(libName);
            return true;
        }
    }
    else
    {
        return false;
    }
}


bool Foam::dlLibraryTable::open
(
    const fileNameList& libNames,
    const bool verbose
)
{
    bool allOpened = !libNames.empty();

    forAll(libNames, i)
    {
        allOpened = open(libNames[i], verbose) && allOpened;
    }

    return allOpened;
}


bool Foam::dlLibraryTable::open
(
    const dictionary& dict,
    const word& libsEntry
)
{
    if (dict.found(libsEntry))
    {
        fileNameList libNames(dict.lookup(libsEntry));

        bool allOpened = !libNames.empty();

        forAll(libNames, i)
        {
            allOpened = open(libNames[i]) && allOpened;
        }

        return allOpened;
    }
    else
    {
        return false;
    }
}


bool Foam::dlLibraryTable::close
(
    const fileName& libName,
    const bool verbose
)
{
    label index = -1;
    forAllReverse(libNames_, i)
    {
        if (libNames_[i] == libName)
        {
            index = i;
            break;
        }
    }

    if (index != -1)
    {
        if (debug)
        {
            InfoInFunction
                << "Closing " << libName
                << " with handle " << uintptr_t(libPtrs_[index]) << endl;
        }

        bool ok = dlClose(libPtrs_[index]);

        libPtrs_[index] = nullptr;
        libNames_[index] = fileName::null;

        if (!ok)
        {
            if (verbose)
            {
                WarningInFunction
                    << "could not close " << libName
                    << endl;
            }

            return false;
        }

        return true;
    }
    return false;
}


void* Foam::dlLibraryTable::findLibrary(const fileName& libName)
{
    label index = -1;
    forAllReverse(libNames_, i)
    {
        if (libNames_[i] == libName)
        {
            index = i;
            break;
        }
    }

    if (index != -1)
    {
        return libPtrs_[index];
    }
    return nullptr;
}


// ************************************************************************* //
