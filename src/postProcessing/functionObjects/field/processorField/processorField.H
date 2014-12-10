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
    Foam::processorField

Group
    grpFieldFunctionObjects

Description
    This function object writes a scalar field whose value is the local
    processor ID.  The output field name is 'processorID'.

    Example of function object specification:
    \verbatim
    processorField1
    {
        type        processorField;
        functionObjectLibs ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: processorField | yes       |
    \endtable

SeeAlso
    Foam::functionObject
    Foam::OutputFilterFunctionObject

SourceFiles
    processorField.C
    IOprocessorField.H

\*---------------------------------------------------------------------------*/

#ifndef processorField_H
#define processorField_H

#include "OFstream.H"
#include "pointFieldFwd.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class processorField Declaration
\*---------------------------------------------------------------------------*/

class processorField
{
protected:

    // Protected data

        //- Name of this set of nearWallFields object
        word name_;

        //- Reference to the database
        const objectRegistry& obr_;

        //- on/off switch
        bool active_;


    // Protected Member Functions

        //- Disallow default bitwise copy construct
        processorField(const processorField&);

        //- Disallow default bitwise assignment
        void operator=(const processorField&);


public:

    //- Runtime type information
    TypeName("processorField");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        processorField
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~processorField();


    // Member Functions

        //- Return name of the processorField object
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the input data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write
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
