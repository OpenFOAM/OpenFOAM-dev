/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Group
    grpUtilitiesFunctionObjects

Description
    This function object provides a general interface to enable dynamic code
    compilation.

    The entries are
        code        : c++; upon functionObject::write()
        codeInclude : include files
        codeOptions : include paths; inserted into EXE_INC in Make/options
        codeLibs    : link line; inserted into LIB_LIBS in Make/options

        codeExecute : c++;upon functionObject::execute();
        codeRead    : c++; upon functionObject::read();
        codeEnd     : c++; upon functionObject::end();
        codeData    : c++; local member data (null constructed);
        codeTimeSet : c++; upon functionObject::timeSet();
        localCode   : c++; local static functions

    Example of function object specification:
    \verbatim
    difference
    {
        functionObjectLibs ("libutilityFunctionObjects.so");

        type coded;
        // Name of on-the-fly generated functionObject
        redirectType writeMagU;
        code
        #{
            // Lookup U
            const volVectorField& U = mesh().lookupObject<volVectorField>("U");
            // Write
            mag(U).write();
        }
    }
    \endverbatim


SeeAlso
    Foam::functionObject
    Foam::OutputFilterFunctionObject
    Foam::codedBase

SourceFiles
    codedFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef codedFunctionObject_H
#define codedFunctionObject_H

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
protected:

    // Protected data

        //- Reference to the time database
        const Time& time_;

        //- Input dictionary
        dictionary dict_;

        word redirectType_;

        string codeData_;
        string codeRead_;
        string codeExecute_;
        string codeEnd_;
        string codeTimeSet_;

        //- Underlying functionObject
        mutable autoPtr<functionObject> redirectFunctionObjectPtr_;


    // Protected Member Functions

        //- get the loaded dynamic libraries
        virtual dlLibraryTable& libs() const;

        //- adapt the context for the current object
        virtual void prepare(dynamicCode &,const dynamicCodeContext&) const;

        // Return a description (type + name) for the output
        virtual string description() const;

        // Clear any redirected objects
        virtual void clearRedirect() const;

        // Get the dictionary to initialize the codeContext
        virtual const dictionary& codeDict() const;


private:

        //- Disallow default bitwise copy construct
        codedFunctionObject(const codedFunctionObject&);

        //- Disallow default bitwise assignment
        void operator=(const codedFunctionObject&);


public:

    //- Runtime type information
    TypeName("coded");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        codedFunctionObject
        (
            const word& name,
            const Time& time,
            const dictionary& dict,
            bool readNow=true   // allow child-classes to avoid compilation
        );


    //- Destructor
    virtual ~codedFunctionObject();


    // Member Functions

        //- Dynamically compiled functionObject
        functionObject& redirectFunctionObject() const;

        //- Called at the start of the time-loop
        virtual bool start();

        //- Called at each ++ or += of the time-loop. forceWrite overrides the
        //  outputControl behaviour.
        virtual bool execute(const bool forceWrite);

        //- Called when Time::run() determines that the time-loop exits.
        //  By default it simply calls execute().
        virtual bool end();

        //- Called when time was set at the end of the Time::operator++
        virtual bool timeSet();

        //- Read and set the function object if its data have changed
        virtual bool read(const dictionary&);

        //- Update mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Move points
        virtual void movePoints(const polyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
