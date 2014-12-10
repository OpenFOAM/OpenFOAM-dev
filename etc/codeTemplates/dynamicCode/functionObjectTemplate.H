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

Description
    Template for use with dynamic code generation of a
    OutputFilter functionObject.

SourceFiles
    functionObjectTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjectTemplate_H
#define functionObjectTemplate_H

#include "stringList.H"
#include "pointField.H"

//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;
class fvMesh;

/*---------------------------------------------------------------------------*\
                         A templated functionObject
\*---------------------------------------------------------------------------*/

class ${typeName}FunctionObject
{
    // Private data

        //- Name of this set of system calls
        word name_;

        //- Registry
        const objectRegistry& obr_;

//{{{ begin codeData
        ${codeData}
//}}} end codeData


    // Private Member Functions

        const objectRegistry& obr() const;

        const fvMesh& mesh() const;

        //- Disallow default bitwise copy construct
        ${typeName}FunctionObject(const ${typeName}FunctionObject&);

        //- Disallow default bitwise assignment
        void operator=(const ${typeName}FunctionObject&);


public:

    //- Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        ${typeName}FunctionObject
        (
            const word& name,
            const objectRegistry& unused,
            const dictionary&,
            const bool loadFromFilesUnused = false
        );


    //- Destructor
    virtual ~${typeName}FunctionObject();


    // Member Functions

        //- Return name of the system call set
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the system calls
        virtual void read(const dictionary&);

        //- Execute the "executeCalls" at each time-step
        virtual void execute();

        //- Execute the "endCalls" at the final time-loop
        virtual void end();

        //- Write, execute the "writeCalls"
        virtual void write();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
