/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Base class for function objects and boundary conditions using dynamic code

SourceFiles
    codedBase.C

\*---------------------------------------------------------------------------*/

#ifndef codedBase_H
#define codedBase_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dynamicCode;
class dynamicCodeContext;
class dlLibraryTable;

/*---------------------------------------------------------------------------*\
                       Class codedBase Declaration
\*---------------------------------------------------------------------------*/

class codedBase
{
    // Private data

        //- Previously loaded library
        mutable fileName oldLibPath_;

    // Private Member Functions

        //- Global loader/unloader function type
        typedef void (*loaderFunctionType)(bool);

        //- Load specified library and execute globalFuncName(true)
        void* loadLibrary
        (
            const fileName& libPath,
            const string& globalFuncName,
            const dictionary& contextDict
        ) const;

        //- Execute globalFuncName(false) and unload specified library
        void unloadLibrary
        (
            const fileName& libPath,
            const string& globalFuncName,
            const dictionary& contextDict
        ) const;

        //- Create library based on the dynamicCodeContext
        void createLibrary(dynamicCode&, const dynamicCodeContext&) const;

        //- Disallow default bitwise copy construct
        codedBase(const codedBase&);

        //- Disallow default bitwise assignment
        void operator=(const codedBase&);


protected:

        //- Update library as required
        void updateLibrary
        (
            const word& redirectType
        ) const;

        //- get the loaded dynamic libraries
        virtual dlLibraryTable& libs() const = 0;

        //- adapt the context for the current object
        virtual void prepare
        (
            dynamicCode&,
            const dynamicCodeContext&
        ) const = 0;

        // Return a description (type + name) for the output
        virtual string description() const = 0;

        // Clear any redirected objects
        virtual void clearRedirect() const = 0;

        // Get the dictionary to initialize the codeContext
        virtual const dictionary& codeDict() const = 0;


public:

    //- Runtime type information
    ClassName("codedBase");


    // Constructors

        //- Construct null
        codedBase();


    //- Destructor
    virtual ~codedBase();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
