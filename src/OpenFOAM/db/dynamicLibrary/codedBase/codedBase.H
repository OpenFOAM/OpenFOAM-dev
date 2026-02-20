/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::codedBase

Description
    Base class for coded functionObjects, fvModels, Function1, Function2

    and boundary conditions

SourceFiles
    codedBase.C

\*---------------------------------------------------------------------------*/

#ifndef codedBase_H
#define codedBase_H

#include "dynamicCode.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class codedBase Declaration
\*---------------------------------------------------------------------------*/

class codedBase
:
    public dynamicCode
{
    // Private Data

        //- Previously loaded library
        mutable fileName oldLibPath_;


    // Private Member Functions

        //- Global loader/unloader function type
        typedef void (*loaderFunctionType)(bool);

        //- Create a code name from the given function name
        static word filterCodeName(const word& name);

        void checkLibrary
        (
            const dictionary& dict,
            const fileName& libPath,
            void* lib,
            const string& uniqueFuncName,
            const bool load
        ) const;

        //- Load specified library and execute uniqueFuncName(true)
        void* loadLibrary
        (
            const dictionary& dict,
            const fileName& libPath,
            const string& uniqueFuncName
        ) const;

        //- Execute uniqueFuncName(false) and unload specified library
        void unloadLibrary
        (
            const dictionary& dict,
            const fileName& libPath,
            const string& uniqueFuncName
        ) const;

        verbatimString expandCodeString
        (
            const dictionary& dict,
            const word& codeKey,
            const word& codeDictVar
        ) const;


public:

    //- Runtime type information
    TypeName("codedBase");


    // Constructors

        //- Construct from name and dictionary
        codedBase
        (
            const word& name,
            const dictionary& dict,
            const wordList& codeKeys,
            const wordList& codeDictVars,
            const word& codeOptionsFileName,
            const wordList& compileFiles,
            const wordList& copyFiles
        );

        //- Construct from dictionary
        codedBase
        (
            const dictionary& dict,
            const wordList& codeKeys,
            const wordList& codeDictVars,
            const word& codeOptionsFileName,
            const wordList& compileFiles,
            const wordList& copyFiles
        );


    //- Destructor
    virtual ~codedBase();


    // Member Functions

        //- Update library from given updated dictionary as required
        //  Returns true if the library was updated, otherwise false
        bool updateLibrary(const dictionary& dict) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const codedBase&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
