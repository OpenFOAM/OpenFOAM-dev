/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::codedFunctionObject

Description
    Provides a general interface to enable dynamic code compilation.

    The entries are:
    \plaintable
       codeInclude | include files
       codeOptions | include paths; inserted into EXE_INC in Make/options
       codeLibs    | link line; inserted into LIB_LIBS in Make/options
       codeData    | c++; local member data (null constructed);
       localCode   | c++; local static functions;
       codeRead    | c++; upon functionObject::read();
       codeFields  | c++; upon functionObject::fields();
       codeExecute | c++; upon functionObject::execute();
       codeWrite   | c++; upon functionObject::write()
       codeEnd     | c++; upon functionObject::end();
    \endplaintable

    Example of function object specification:
    \verbatim
    writeMagU
    {
        libs        ("libutilityFunctionObjects.so");

        type coded;

        codeWrite
        #{
            // Lookup U
            const volVectorField& U = mesh().lookupObject<volVectorField>("U");
            // Write
            mag(U)().write();
        #};
    }
    \endverbatim

See also
    Foam::functionObject
    Foam::codedBase

SourceFiles
    codedFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_codedFunctionObject_H
#define functionObjects_codedFunctionObject_H

#include "functionObject.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class codedFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class codedFunctionObject
:
    public functionObject,
    public codedBase
{
    // Private Data

        //- Reference to the time database
        const Time& time_;

        //- The dynamically generated functionObject pointer
        mutable autoPtr<functionObject> redirectFunctionObjectPtr_;


    // Private Member Functions

        //- Get the keywords associated with source code
        virtual wordList codeKeys() const;

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;

        //- Clear any redirected objects
        virtual void clearRedirect() const;

        //- Compile, link and return the now coded functionObject
        functionObject& redirectFunctionObject() const;


public:

    //- Runtime type information
    TypeName("coded");


    // Constructors

        //- Construct from Time and dictionary
        codedFunctionObject
        (
            const word& name,
            const Time& time,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        codedFunctionObject(const codedFunctionObject&) = delete;


    //- Destructor
    virtual ~codedFunctionObject();


    // Member Functions

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Called at each ++ or += of the time-loop.
        //  foamPostProcess overrides the usual executeControl behaviour and
        //  forces execution (used in post-processing mode)
        virtual bool execute();

        //- Called at each ++ or += of the time-loop.
        //  foamPostProcess overrides the usual writeControl behaviour and
        //  forces writing always (used in post-processing mode)
        virtual bool write();

        //- Called when Time::run() determines that the time-loop exits.
        //  By default it simply calls execute().
        virtual bool end();

        //- Read and set the function object if its data have changed
        virtual bool read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const codedFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
