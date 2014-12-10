/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Foam::fieldValues::fieldValueDelta

Group
    grpFieldFunctionObjects

Description
    This function object provides a differencing option between two 'field
    value' function objects.

    Example of function object specification:
    \verbatim
    fieldValueDelta1
    {
        type            fieldValueDelta;
        functionObjectLibs ("libfieldFunctionObjects.so");
        operation       subtract;

        fieldValue1
        {
            ...
        }
        fieldValue2
        {
            ...
        }
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: fieldValueDelta   | yes    |
    \endtable

    \linebreak
    The \c operation is one of:
    \plaintable
       add           | add
       subtract      | subtract
       min           | minimum
       max           | maximum
       average       | average
    \endplaintable
SeeAlso
    Foam::fieldValue

SourceFiles
    fieldValueDelta.C

\*---------------------------------------------------------------------------*/

#ifndef fieldValueDelta_H
#define fieldValueDelta_H

#include "functionObjectFile.H"
#include "fieldValue.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fieldValues
{

/*---------------------------------------------------------------------------*\
                       Class fieldValueDelta Declaration
\*---------------------------------------------------------------------------*/

class fieldValueDelta
:
    public functionObjectFile
{
public:
        //- Operation type enumeration
        enum operationType
        {
            opAdd,
            opSubtract,
            opMin,
            opMax,
            opAverage
        };

        //- Operation type names
        static const NamedEnum<operationType, 5> operationTypeNames_;


private:

    // Private data

        //- Name of this fieldValue object
        word name_;

        //- Database this class is registered to
        const objectRegistry& obr_;

        //- Flag to indicate to load from files
        bool loadFromFiles_;

        //- Switch to send output to Info as well as to file
        Switch log_;

        //- Operation to apply to values
        operationType operation_;

        //- Field value source object 1
        autoPtr<fieldValue> source1Ptr_;

        //- Field value source object 2
        autoPtr<fieldValue> source2Ptr_;


    // Private Member Functions

        //- Templated function to process common fields
        template<class Type>
        void processFields(bool& found);

        //- Templated function to apply the operation
        template<class Type>
        Type applyOperation(const Type& value1, const Type& value2) const;


protected:

    // Functions to be over-ridden from IOoutputFilter class

        //- Output file header information
        virtual void writeFileHeader(const label i);


public:

    //- Run-time type information
    TypeName("fieldValueDelta");


    //- Construct from components
    fieldValueDelta
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        const bool loadFromFiles = false
    );


    //- Destructor
    virtual ~fieldValueDelta();


    // Public Member Functions

        // Function object functions

            //- Read from dictionary
            virtual void read(const dictionary&);

            //- Calculate and write
            virtual void write();

            //- Execute
            virtual void execute();

            //- Execute the at the final time-loop, currently does nothing
            virtual void end();

             //- Called when time was set at the end of the Time::operator++
            virtual void timeSet();

            //- Update mesh
            virtual void updateMesh(const mapPolyMesh&);

            //- Move points
            virtual void movePoints(const polyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldValues
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldValueDeltaTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
