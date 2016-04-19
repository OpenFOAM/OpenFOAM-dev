/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::FUNCTIONOBJECT

Group

Description
    This function object...

    Example of function object specification:
    \verbatim
    FUNCTIONOBJECT1
    {
        type           FUNCTIONOBJECT;
        functionObjectLibs ("libFUNCTIONOBJECTFunctionObject.so");
        ...
        wordData       someWord;
        scalarData     1.0;
        labelData      1;
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description               | Required | Default value
        type         | type name: FUNCTIONOBJECT | yes      |
        wordData     | some word option...       | no       | defaultWord
        scalarData   | some scalar value...      | yes      |
        labelData    | some label value...       | yes      |
    \endtable

SourceFiles
    FUNCTIONOBJECT.C
    IOFUNCTIONOBJECT.H

\*---------------------------------------------------------------------------*/

#ifndef FUNCTIONOBJECT_H
#define FUNCTIONOBJECT_H

#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                   Class FUNCTIONOBJECT Declaration
\*---------------------------------------------------------------------------*/

class FUNCTIONOBJECT
{
    // Private data

        //- Name of this set of FUNCTIONOBJECT
        word name_;

        //- Refererence to Db
        const objectRegistry& obr_;

        // Read from dictionary

            //- word
            word wordData_;

            //- scalar
            scalar scalarData_;

            //- label
            label labelData_;


    // Private Member Functions


        //- Disallow default bitwise copy construct
        FUNCTIONOBJECT(const FUNCTIONOBJECT&);

        //- Disallow default bitwise assignment
        void operator=(const FUNCTIONOBJECT&);


public:

    //- Runtime type information
    TypeName("FUNCTIONOBJECT");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        FUNCTIONOBJECT
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~FUNCTIONOBJECT();


    // Member Functions

        //- Return name of the FUNCTIONOBJECT
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the FUNCTIONOBJECT data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write the FUNCTIONOBJECT
        virtual void write();

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
