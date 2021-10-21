/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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
    Foam::functionObjects::fieldValues::fieldValueDelta

Description
    Provides a differencing option between two 'field value' function objects.

    Example of function object specification:
    \verbatim
    fieldValueDelta1
    {
        type            fieldValueDelta;
        libs            ("libfieldFunctionObjects.so");
        operation       subtract;

        region1
        {
            ...
        }

        region2
        {
            ...
        }
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: fieldValueDelta   | yes    |
    \endtable

    The \c operation is one of:
    \plaintable
       add           | add
       subtract      | subtract
       min           | minimum
       max           | maximum
       average       | average
    \endplaintable

See also
    Foam::functionObject
    Foam::functionObject::fieldValue
    Foam::functionObjects::regionFunctionObject
    Foam::functionObjects::logFiles

SourceFiles
    fieldValueDelta.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldValueDelta_H
#define functionObjects_fieldValueDelta_H

#include "regionFunctionObject.H"
#include "logFiles.H"
#include "fieldValue.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{

/*---------------------------------------------------------------------------*\
                       Class fieldValueDelta Declaration
\*---------------------------------------------------------------------------*/

class fieldValueDelta
:
    public regionFunctionObject,
    public logFiles
{
public:
        //- Operation type enumeration
        enum class operationType
        {
            add,
            subtract,
            min,
            max,
            average
        };

        //- Operation type names
        static const NamedEnum<operationType, 5> operationTypeNames_;


private:

    // Private Data

        //- Operation to apply to values
        operationType operation_;

        //- Field value region object 1
        autoPtr<fieldValue> region1Ptr_;

        //- Field value region object 2
        autoPtr<fieldValue> region2Ptr_;


    // Private Member Functions

        //- Templated function to process common fields
        template<class Type>
        void processFields(bool& found);

        //- Templated function to apply the operation
        template<class Type>
        Type applyOperation(const Type& value1, const Type& value2) const;


protected:

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(const label i);


public:

    //- Run-time type information
    TypeName("fieldValueDelta");


    // Constructors

        //- Construct from Time and dictionary
        fieldValueDelta
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fieldValueDelta();


    // Public Member Functions

        //- Read from dictionary
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Do nothing
        virtual bool execute();

        //- Calculate and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldValues
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldValueDeltaTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
