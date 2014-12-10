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
    Foam::fieldMinMax

Group
    grpFieldFunctionObjects

Description
    This function object calculates the value and position of scalar minimim
    and maximum for a list of user-specified fields.  For variables with a rank
    greater than zero, either the min/max of a component value or the magnitude
    is reported.  When operating in parallel, the processor owning the value
    is also given.

    Example of function object specification:
    \verbatim
    fieldMinMax1
    {
        type        fieldMinMax;
        functionObjectLibs ("libfieldFunctionObjects.so");
        ...
        write       yes;
        log         yes;
        mode        magnitude;
        fields
        (
            U
            p
        );
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: fieldMinMax  | yes         |
        write        | write min/max data to file |  no      | yes
        log          | write min/max data to standard output | no | no
        mode         | calculation mode: magnitude or component | no | magnitude
    \endtable

    Output data is written to the file \<timeDir\>/fieldMinMax.dat

SeeAlso
    Foam::functionObject
    Foam::OutputFilterFunctionObject

SourceFiles
    fieldMinMax.C
    IOfieldMinMax.H

\*---------------------------------------------------------------------------*/

#ifndef fieldMinMax_H
#define fieldMinMax_H

#include "functionObjectFile.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "OFstream.H"
#include "Switch.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class fieldMinMax Declaration
\*---------------------------------------------------------------------------*/

class fieldMinMax
:
    public functionObjectFile
{
public:

    enum modeType
    {
        mdMag,
        mdCmpt
    };

protected:

    // Protected data

        //- Mode type names
        static const NamedEnum<modeType, 2> modeTypeNames_;

        //- Name of this set of field min/max
        //  Also used as the name of the output directory
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Switch to send output to Info as well
        Switch log_;

        //- Mode for min/max - only applicable for ranks > 0
        modeType mode_;

        //- Fields to assess min/max
        wordList fieldSet_;


    // Private Member Functions

        //- Helper function to write the output
        template<class Type>
        void output
        (
            const word& fieldName,
            const word& outputName,
            const vector& minC,
            const vector& maxC,
            const label minProcI,
            const label maxProcI,
            const Type& minValue,
            const Type& maxValue
        );

        //- Disallow default bitwise copy construct
        fieldMinMax(const fieldMinMax&);

        //- Disallow default bitwise assignment
        void operator=(const fieldMinMax&);

        //- Output file header information
        virtual void writeFileHeader(const label i);


public:

    //- Runtime type information
    TypeName("fieldMinMax");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        fieldMinMax
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~fieldMinMax();


    // Member Functions

        //- Return name of the set of field min/max
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the field min/max data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the field min/max
        template<class Type>
        void calcMinMaxFields
        (
            const word& fieldName,
            const modeType& mode
        );

        //- Write the fieldMinMax
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

#ifdef NoRepository
#   include "fieldMinMaxTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
