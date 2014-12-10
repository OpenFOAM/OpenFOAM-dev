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

\*---------------------------------------------------------------------------*/

#include "dlLibraryTable.H"
#include "dictionary.H"
#include "fileNameList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TablePtr>
bool Foam::dlLibraryTable::open
(
    const dictionary& dict,
    const word& libsEntry,
    const TablePtr& tablePtr
)
{
    if (dict.found(libsEntry))
    {
        fileNameList libNames(dict.lookup(libsEntry));

        bool allOpened = (libNames.size() > 0);

        forAll(libNames, i)
        {
            const fileName& libName = libNames[i];

            label nEntries = 0;

            if (tablePtr)
            {
                nEntries = tablePtr->size();
            }

            bool opened = dlLibraryTable::open(libName);
            allOpened = opened && allOpened;

            if (!opened)
            {
                WarningIn
                (
                    "dlLibraryTable::open"
                    "(const dictionary&, const word&, "
                    "const TablePtr&)"
                )   << "Could not open library " << libName
                    << endl << endl;
            }
            else if (debug && (!tablePtr || tablePtr->size() <= nEntries))
            {
                WarningIn
                (
                    "dlLibraryTable::open"
                    "(const dictionary&, const word&, "
                    "const TablePtr&)"
                )   << "library " << libName
                    << " did not introduce any new entries"
                    << endl << endl;
            }
        }

        return allOpened;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
