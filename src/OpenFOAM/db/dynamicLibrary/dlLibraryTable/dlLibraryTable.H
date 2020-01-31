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

Class
    Foam::dlLibraryTable

Description
    A table of dynamically loaded libraries

SourceFiles
    dlLibraryTable.C

\*---------------------------------------------------------------------------*/

#ifndef dlLibraryTable_H
#define dlLibraryTable_H

#include "label.H"
#include "DynamicList.H"
#include "fileNameList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class dlLibraryTable Declaration
\*---------------------------------------------------------------------------*/

class dlLibraryTable
{
    // Private Member Functions

        //- List of pointers to the loaded dynamic libraries
        DynamicList<void*> libPtrs_;

        //- List of names of the loaded dynamic libraries
        DynamicList<fileName> libNames_;


public:

    // Declare name of the class and its debug switch
    ClassName("dlLibraryTable");


    // Constructors

        //- Construct null
        dlLibraryTable();

        //- Construct from dictionary and name of 'libs' entry giving
        //  the libraries to load
        dlLibraryTable(const dictionary&, const word&);

        //- Disallow default bitwise copy construction
        dlLibraryTable(const dlLibraryTable&) = delete;


    //- Destructor
    ~dlLibraryTable();


    // Member Functions

        //- Open the named library, optionally with warnings if problems occur
        bool open(const fileName& libName, const bool verbose = true);

        //- Open the named libraries, optionally with warnings if problems occur
        bool open(const fileNameList& libNames, const bool verbose = true);

        //- Open all the libraries listed in the 'libsEntry' entry in the
        //  given dictionary if present
        bool open(const dictionary&, const word& libsEntry);

        //- Open all the libraries listed in the 'libsEntry' entry in the
        //  given dictionary if present and check the additions
        //  to the given constructor table
        template<class TablePtr>
        bool open
        (
            const dictionary&,
            const word& libsEntry,
            const TablePtr& tablePtr
        );

        //- Find the handle of the named library
        void* findLibrary(const fileName& libName);

        //- Close the named library, optionally with warnings if problems occur
        bool close(const fileName& name, const bool verbose = true);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dlLibraryTable&) = delete;
};


//- Table of loaded dynamic libraries
extern dlLibraryTable libs;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dlLibraryTableTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
