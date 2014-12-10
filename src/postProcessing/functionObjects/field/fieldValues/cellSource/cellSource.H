/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::fieldValues::cellSource

Group
    grpFieldFunctionObjects

Description
    This function object provides a 'cell source' variant of the fieldValues
    function object.  Given a list of user-specified fields and a selection
    of mesh cells, a number of operations can be performed, such as sums,
    averages and integrations.


    Example of function object specification:
    \verbatim
    cellSource1
    {
        type            cellSource;
        functionObjectLibs ("libfieldFunctionObjects.so");
        ...
        log             true;
        valueOutput     true;
        source          cellZone;
        sourceName      c0;
        operation       volAverage;
        weightField     alpha1;
        fields
        (
            p
            U
        );
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: cellSource   | yes         |
        log          | write data to standard output | no    | no
        valueOutput  | write the raw output values | yes     |
        source       | cell source: see below  | yes         |
        sourceName   | name of cell source if required  | no |
        operation    | operation to perform    | yes         |
        weightField  | name of field to apply weighting | no |
        fields       | list of fields to operate on | yes    |
    \endtable

    \linebreak
    Where \c source is defined by
    \plaintable
        cellZone     | requires a 'sourceName' entry to specify the cellZone
        all          | all cells
    \endplaintable

    \linebreak
    The \c operation is one of:
    \plaintable
       none          | no operation
       sum           | sum
       average       | ensemble average
       weightedAverage | weighted average
       volAverage    | volume weighted average
       volIntegrate  | volume integral
       min           | minimum
       max           | maximum
       CoV           | coefficient of variation: standard deviation/mean
    \endplaintable

SeeAlso
    Foam::fieldValues
    Foam::functionObject
    Foam::OutputFilterFunctionObject

SourceFiles
    cellSource.C

\*---------------------------------------------------------------------------*/

#ifndef cellSource_H
#define cellSource_H

#include "NamedEnum.H"
#include "fieldValue.H"
#include "labelList.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldValues
{

/*---------------------------------------------------------------------------*\
                         Class cellSource Declaration
\*---------------------------------------------------------------------------*/

class cellSource
:
    public fieldValue
{

public:

    // Public data types

        //- Source type enumeration
        enum sourceType
        {
            stCellZone,
            stAll
        };

        //- Source type names
        static const NamedEnum<sourceType, 2> sourceTypeNames_;


        //- Operation type enumeration
        enum operationType
        {
            opNone,
            opSum,
            opAverage,
            opWeightedAverage,
            opVolAverage,
            opVolIntegrate,
            opMin,
            opMax,
            opCoV
        };

        //- Operation type names
        static const NamedEnum<operationType, 9> operationTypeNames_;


private:

    // Private Member Functions

        //- Set cells to evaluate based on a cell zone
        void setCellZoneCells();

        //- Set cells to evaluate based on a patch
        void setPatchCells();


protected:

    // Protected data

        //- Source type
        sourceType source_;

        //- Operation to apply to values
        operationType operation_;

        //- Global number of cells
        label nCells_;

        //- Local list of cell IDs
        labelList cellId_;

        //- Weight field name - only used for opWeightedAverage mode
        word weightFieldName_;


    // Protected Member Functions

        //- Initialise, e.g. cell addressing
        void initialise(const dictionary& dict);

        //- Return true if the field name is valid
        template<class Type>
        bool validField(const word& fieldName) const;

        //- Insert field values into values list
        template<class Type>
        tmp<Field<Type> > setFieldValues
        (
            const word& fieldName,
            const bool mustGet = false
        ) const;

        //- Apply the 'operation' to the values
        template<class Type>
        Type processValues
        (
            const Field<Type>& values,
            const scalarField& V,
            const scalarField& weightField
        ) const;

        //- Output file header information
        virtual void writeFileHeader(const label i);


public:

    //- Run-time type information
    TypeName("cellSource");


    //- Construct from components
    cellSource
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        const bool loadFromFiles = false
    );


    //- Destructor
    virtual ~cellSource();


    // Public Member Functions

        // Access

            //- Return the source type
            inline const sourceType& source() const;

            //- Return the local list of cell IDs
            inline const labelList& cellId() const;


        // Function object functions

            //- Read from dictionary
            virtual void read(const dictionary&);

            //- Calculate and write
            virtual void write();

            //- Templated helper function to output field values
            template<class Type>
            bool writeValues(const word& fieldName);

            //- Filter a field according to cellIds
            template<class Type>
            tmp<Field<Type> > filterField(const Field<Type>& field) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldValues
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellSourceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cellSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
