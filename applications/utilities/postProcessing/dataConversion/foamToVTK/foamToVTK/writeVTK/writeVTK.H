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
    Foam::functionObjects::writeVTK

Group
    grpUtilitiesFunctionObjects

Description
    This functionObject writes objects registered to the database in VTK format
    using the foamToVTK library.

    Currently only the writing of the cell-values of volFields is supported but
    support for other field types, patch fields, Lagrangian data etc. will be
    added.

    Example of function object specification:
    \verbatim
        writeVTK1
        {
            type        writeVTK;
            functionObjectLibs ("libIOFunctionObjects.so");
            ...
            objectNames (obj1 obj2);
        }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: writeVTK     | yes         |
        objectNames  | objects to write        | yes         |
    \endtable

SeeAlso
    Foam::functionObject
    Foam::OutputFilterFunctionObject

SourceFiles
    writeVTK.C
    IOwriteVTK.H

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeVTK_H
#define functionObjects_writeVTK_H

#include "wordReList.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class writeVTK Declaration
\*---------------------------------------------------------------------------*/

class writeVTK
{
    // Private data

        //- Name of this set of writeVTK
        word name_;

        //- Refererence to Db
        const objectRegistry& obr_;

        //- Names of objects
        wordReList objectNames_;


    // Private Member Functions

        template<class GeoField>
        UPtrList<const GeoField> lookupFields() const;

        //- Disallow default bitwise copy construct
        writeVTK(const writeVTK&);

        //- Disallow default bitwise assignment
        void operator=(const writeVTK&);


public:

    //- Runtime type information
    TypeName("writeVTK");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        writeVTK
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~writeVTK();


    // Member Functions

        //- Return name of the writeVTK
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the writeVTK data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write the writeVTK
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "writeVTKTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
