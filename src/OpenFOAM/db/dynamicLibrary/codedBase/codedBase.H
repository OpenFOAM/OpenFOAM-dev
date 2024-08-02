/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "dynamicCodeContext.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dynamicCode;

/*---------------------------------------------------------------------------*\
                          Class codedBase Declaration
\*---------------------------------------------------------------------------*/

class codedBase
{
    // Private Member Data

        //- Name of the dynamically generated CodedType
        const word codeName_;

        dynamicCodeContext codeContext_;

        //- Previously loaded library
        mutable fileName oldLibPath_;


    // Private Member Functions

        //- Global loader/unloader function type
        typedef void (*loaderFunctionType)(bool);

        //- Create a code name from the given function name
        static word codeName(const word& name);

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

        verbatimString expandCodeString
        (
            const word& codeKey,
            const word& codeDictVar,
            const dictionary& dict
        ) const;

        //- Create library based on the dynamicCodeContext
        void createLibrary
        (
            const dictionary& dict,
            dynamicCode&,
            const dynamicCodeContext&
        ) const;

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const = 0;


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
            const wordList& codeDictVars
        );

        //- Construct from dictionary
        codedBase
        (
            const dictionary& dict,
            const wordList& codeKeys,
            const wordList& codeDictVars
        );

        //- Copy constructor
        codedBase(const codedBase& cb);


    //- Destructor
    virtual ~codedBase();


    // Member Functions

        //- Name of the dynamically generated CodedType
        const word& codeName() const;

        //- Return a description (type + name) for the output
        string description() const;

        word codeTemplateC(const word& baseTypeName) const;

        word codeTemplateH(const word& baseTypeName) const;

        //- Update library from given updated dictionary as required
        //  Returns true if the library was updated, otherwise false
        bool updateLibrary(const dictionary& dict) const;

        //- Read the dictionary and update the code
        void read(const dictionary& dict);

        //- Write the code for restart
        void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const codedBase&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
