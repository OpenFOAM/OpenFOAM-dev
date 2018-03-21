/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::fileOperationInitialise

\*---------------------------------------------------------------------------*/

#ifndef fileOperationInitialise_H
#define fileOperationInitialise_H

#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileOperations
{

/*---------------------------------------------------------------------------*\
                   Class fileOperationInitialise Declaration
\*---------------------------------------------------------------------------*/

class fileOperationInitialise
{
public:

        //- Runtime type information
        TypeName("fileOperationInitialise");


    // Constructors

        //- Construct components
        fileOperationInitialise(int& argc, char**& argv);


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fileOperationInitialise,
            word,
            (
                int& argc, char**& argv
            ),
            (argc, argv)
        );


    // Selectors

        //- Select type
        static autoPtr<fileOperationInitialise> New
        (
            const word& type, int& argc, char**& argv
        );


    //- Destructor
    virtual ~fileOperationInitialise();


    // Member Functions

        //- Needs threading
        virtual bool needsThreading() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileOperations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
