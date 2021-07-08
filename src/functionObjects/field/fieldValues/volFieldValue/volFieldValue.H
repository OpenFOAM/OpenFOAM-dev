/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::functionObjects::fieldValues::volFieldValue

Description
    Provides a 'volRegion' specialisation of the fieldValue function object.

    Given a list of user-specified fields and a 'volRegion', a number of
    operations can be performed, such as sums, averages and integrations.

    Example of function object specification:
    \verbatim
    volFieldValue1
    {
        type            volFieldValue;
        libs            ("libfieldFunctionObjects.so");

        log             true;
        writeControl    writeTime;
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
        Property     | Description              | Required    | Default value
        type         | Type name: volFieldValue            | yes   |
        log          | Write data to standard output       | no    | no
        writeFields  | Write the region field values       | yes   |
        writeLocation| Write the location (if available)   | no    | no
        regionType   | volRegion type: see below           | yes   |
        name         | Name of volRegion if required       | no    |
        operation    | Operation to perform                | yes   |
        weightField  | Name of field to apply weighting    | no    |
        weightFields | Names of fields to apply weighting  | no    |
        fields       | List of fields to operate on        | yes   |
    \endtable

    Where \c regionType is defined by
    \plaintable
        cellZone     | requires a 'name' entry to specify the cellZone
        all          | all cells
    \endplaintable

    The \c operation is one of:
    \plaintable
       none                 | No operation
       sum                  | Sum
       sumMag               | Sum of component magnitudes
       average              | Ensemble average
       volAverage           | Volume weighted average
       volIntegrate         | Volume integral
       min                  | Minimum
       max                  | Maximum
       minMag               | Minimum magnitude
       maxMag               | Maximum magnitude
       CoV                  | Coefficient of variation: standard deviation/mean
    \endplaintable

See also
    Foam::functionObjects::fieldValues::fieldValue
    Foam::functionObjects::volRegion
    Foam::functionObject

SourceFiles
    volFieldValue.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_volFieldValue_H
#define functionObjects_volFieldValue_H

#include "fieldValue.H"
#include "volRegion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{

/*---------------------------------------------------------------------------*\
                        Class volFieldValue Declaration
\*---------------------------------------------------------------------------*/

class volFieldValue
:
    public fieldValue,
    public volRegion
{

public:

    // Public data types

        //- Operation type enumeration
        enum class operationType
        {
            none,
            sum,
            sumMag,
            average,
            volAverage,
            volIntegrate,
            min,
            max,
            minMag,
            maxMag,
            CoV
        };

        //- Operation type names
        static const NamedEnum<operationType, 11> operationTypeNames_;


    // Public classes

        //- Forward declare the result structure
        template<class Type>
        struct Result;


protected:

    // Protected data

        //- Write the location if available for this operation - optional
        Switch writeLocation_;

        //- Operation to apply to values
        operationType operation_;

        //- Weight field names
        wordList weightFieldNames_;

        //- Scale factor - optional
        scalar scaleFactor_;


    // Protected Member Functions

        //- Initialise, e.g. cell addressing
        void initialise(const dictionary& dict);

        //- Return true if the field name is valid
        template<class Type>
        bool validField(const word& fieldName) const;

        //- Insert field values into values list
        template<class Type>
        tmp<Field<Type>> getFieldValues(const word& fieldName) const;

        //- Apply a comparison operation (min/max) to the field magnitude,
        //  returning the limiting value, its index and processor index
        template<class Type, class Op>
        void opMag
        (
            const Field<Type>& values,
            Result<scalar>& result,
            const Op& op
        ) const;

        //- Apply the operation to the values, and return true if successful.
        //  Does nothing unless overloaded below.
        template<class Type, class ResultType>
        bool processValues
        (
            const Field<Type>& values,
            const scalarField& weights,
            const scalarField& V,
            Result<ResultType>& result
        ) const;

        //- Apply Type -> Type operation to the values. Calls
        //  processValuesTypeType.
        template<class Type>
        bool processValues
        (
            const Field<Type>& values,
            const scalarField& weights,
            const scalarField& V,
            Result<Type>& result
        ) const;

        //- Apply Type -> scalar operation to the values
        template<class Type>
        bool processValues
        (
            const Field<Type>& values,
            const scalarField& weights,
            const scalarField& V,
            Result<scalar>& result
        ) const;

        //- Apply scalar -> scalar operation to the values. Calls
        //  processValuesTypeType.
        bool processValues
        (
            const Field<scalar>& values,
            const scalarField& weights,
            const scalarField& V,
            Result<scalar>& result
        ) const;

        //- Apply a Type -> Type operation to the values
        template<class Type>
        bool processValuesTypeType
        (
            const Field<Type>& values,
            const scalarField& weights,
            const scalarField& V,
            Result<Type>& result
        ) const;

        //- Output file header location information for a given type
        template<class Type>
        void writeFileHeaderLocation();

        //- Output file header information
        virtual void writeFileHeader(const label i);


public:

    //- Run-time type information
    TypeName("volFieldValue");


    // Constructors

        //- Construct from name, Time and dictionary
        volFieldValue
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from name, objectRegistry and dictionary
        volFieldValue
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~volFieldValue();


    // Public Member Functions

        //- Templated helper function to output field values
        template<class Type>
        bool writeValues
        (
            const word& fieldName,
            const scalarField& weights,
            const scalarField& V
        );

        //- Templated helper function to output field values
        template<class Type, class ResultType>
        bool writeValues
        (
            const word& fieldName,
            const Field<Type>& values,
            const scalarField& weights,
            const scalarField& V
        );

        //- Filter a field according to cellIds
        template<class Type>
        tmp<Field<Type>> filterField(const Field<Type>& field) const;

        //- Read from dictionary
        virtual bool read(const dictionary&);

        //- Calculate and write
        virtual bool write();
};


template<>
void volFieldValue::writeFileHeaderLocation<scalar>();


/*---------------------------------------------------------------------------*\
                      Class volFieldValue::Result Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
struct volFieldValue::Result
{
    Type value;
    label celli;
    label proci;
    point cc;
};


template<class Type>
inline Istream& operator>>
(
    Istream& is,
    volFieldValue::Result<Type>& result
)
{
    return is >> result.value >> result.celli >> result.proci >> result.cc;
}


template<class Type>
inline Ostream& operator<<
(
    Ostream& os,
    const volFieldValue::Result<Type>& result
)
{
    return os << result.value << result.celli << result.proci << result.cc;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldValues
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "volFieldValueTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
