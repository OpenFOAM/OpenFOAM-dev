/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::functionObjects::fieldValues::volRegion

Group
    grpFieldFunctionObjects

Description
    This function object provides a 'cell region' variant of the fieldValues
    function object.  Given a list of user-specified fields and a selection
    of mesh cells, a number of operations can be performed, such as sums,
    averages and integrations.


    Example of function object specification:
    \verbatim
    volRegion1
    {
        type            volRegion;
        libs ("libfieldFunctionObjects.so");
        ...
        log             true;
        writeFields     true;
        regionType      cellZone;
        name            c0;
        operation       volAverage;
        weightField     alpha1;
        fields
        (
            p
            U
        );
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: volRegion   | yes         |
        log          | Write data to standard output | no    | no
        writeFields  | Write the region field values  | yes     |
        writeVolume  | Write the volume of the volRegion | no |
        regionType   | cell regionType: see below  | yes         |
        name         | name of cell regionType if required  | no |
        operation    | operation to perform    | yes         |
        weightField  | name of field to apply weighting | no |
        fields       | list of fields to operate on | yes    |
    \endtable

    Where \c regionType is defined by
    \plaintable
        cellZone     | requires a 'name' entry to specify the cellZone
        all          | all cells
    \endplaintable

    The \c operation is one of:
    \plaintable
       none          | no operation
       sum           | sum
       sumMag        | sum of component magnitudes
       average       | ensemble average
       weightedAverage | weighted average
       volAverage    | volume weighted average
       weightedVolAverage | weighted volume average
       volIntegrate  | volume integral
       min           | minimum
       max           | maximum
       CoV           | coefficient of variation: standard deviation/mean
    \endplaintable

See also
    Foam::fieldValues
    Foam::functionObject

SourceFiles
    volRegion.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_volRegion_H
#define functionObjects_volRegion_H

#include "fieldValue.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{

/*---------------------------------------------------------------------------*\
                         Class volRegion Declaration
\*---------------------------------------------------------------------------*/

class volRegion
:
    public fieldValue
{

public:

    // Public data types

        //- region type enumeration
        enum regionTypes
        {
            stCellZone,
            stAll
        };

        //- region type names
        static const NamedEnum<regionTypes, 2> regionTypeNames_;


        //- Operation type enumeration
        enum operationType
        {
            opNone,
            opSum,
            opSumMag,
            opAverage,
            opWeightedAverage,
            opVolAverage,
            opWeightedVolAverage,
            opVolIntegrate,
            opMin,
            opMax,
            opCoV
        };

        //- Operation type names
        static const NamedEnum<operationType, 11> operationTypeNames_;


private:

    // Private Member Functions

        //- Set cells to evaluate based on a cell zone
        void setCellZoneCells();

        //- Set cells to evaluate based on a patch
        void setPatchCells();

        //- Calculate and return volume of the volRegion: sum(V)
        scalar volume() const;


protected:

    // Protected data

        //- region type
        regionTypes regionType_;

        //- Operation to apply to values
        operationType operation_;

        //- Global number of cells
        label nCells_;

        //- Local list of cell IDs
        labelList cellId_;

        //- Weight field name - only used for opWeightedAverage mode
        word weightFieldName_;

        //- Volume of the volRegion
        scalar volume_;

        //- Optionally write the volume of the volRegion
        bool writeVolume_;


    // Protected Member Functions

        //- Initialise, e.g. cell addressing
        void initialise(const dictionary& dict);

        //- Return true if the field name is valid
        template<class Type>
        bool validField(const word& fieldName) const;

        //- Insert field values into values list
        template<class Type>
        tmp<Field<Type>> setFieldValues
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
    TypeName("volRegion");


    // Constructors

        //- Construct from name, Time and dictionary
        volRegion
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from name, objectRegistry and dictionary
        volRegion
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~volRegion();


    // Public Member Functions

        //- Return the region type
        inline const regionTypes& regionType() const;

        //- Return the local list of cell IDs
        inline const labelList& cellId() const;

        //- Templated helper function to output field values
        template<class Type>
        bool writeValues(const word& fieldName);

        //- Filter a field according to cellIds
        template<class Type>
        tmp<Field<Type>> filterField(const Field<Type>& field) const;

        //- Read from dictionary
        virtual bool read(const dictionary&);

        //- Calculate and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldValues
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volRegionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "volRegionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
