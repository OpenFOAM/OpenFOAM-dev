/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::functionObjects::LagrangianFieldValue

Description
    Function to log a single reduced quantity generated from the values in a
    Lagrangian field; e.g., sums, averages, maximums and minimums.

Usage
    \table
        Property      | Description                  | Required? | Default
        Lagrangian    | Name of the Lagrangian mesh  | yes       |
        field         | Field to operate on          | if fields not specified |
        fields        | List of fields to operate on | if field not specified |
        weightField   | Field with which to weight the distribution | no | none
        weightFields  | List of fields with which to \
                        weight the distribution      | no        | none
        operation     | The operation with which to \
                        combine Lagrangian values    | yes       |
        writeLocation | Whether or not to write the location | no | false
    \endtable

    Where \c operation is one of
    \plaintable
        sum     | Sum
        average | Ensemble average
        min     | Minimum (component minimum if a higher rank type)
        max     | Maximum (component maximum if a higher rank type)
        minMag  | Minimum magnitude
        maxMag  | Maximum magnitude
    \endplaintable

    Example specification to generate the total mass:
    \verbatim
    LagrangianFieldValue1
    {
        type            LagrangianFieldValue;
        libs            ("libLagrangianFunctionObjects.so");
        Lagrangian      cloud;
        field           m;
        weightField     number;
    }
    \endverbatim

SourceFiles
    LagrangianFieldValue.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianFieldValue_H
#define LagrangianFieldValue_H

#include "LagrangianMeshFunctionObject.H"
#include "logFiles.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class LagrangianFieldValue Declaration
\*---------------------------------------------------------------------------*/

class LagrangianFieldValue
:
    public LagrangianMeshFunctionObject,
    public logFiles
{
public:

    // Public Data Types

        //- Operation type enumeration
        enum class operationType
        {
            sum,
            average,
            min,
            max,
            minMag,
            maxMag
        };

        //- Operation type names
        static const NamedEnum<operationType, 6> operationTypeNames_;


private:

    // Private Data

        //- List of fields
        wordList fields_;

        //- List of weight fields
        wordList weightFields_;

        //- Operation to apply to values
        operationType operation_;

        //- Write the location if available for this operation
        Switch writeLocation_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Write a name
        template<class Type>
        void writeName(const word& name);

        //- Write a location name
        template<class Type, class LocationType>
        void writeLocationName(const word& name, const word& locationName);

        //- Write a name and a location name if we are writing the location
        template<class Type>
        void writeNameAndLocationNames(const word& name);

        //- Write a value
        template<class Type>
        void writeValue(const Type& value);

        //- Write a location value
        template<class Type, class LocationType>
        void writeLocationValue
        (
            const FixedList<LocationType, pTraits<Type>::nComponents>& value
        );

        //- Write a value and location values if we are writing the location
        template<class Type, class Op>
        void writeValueAndLocationValues
        (
            const tmp<Field<Type>>& tField,
            const scalar emptyValue,
            const Op& op
        );

        //- Multiply the argument by the given weight field. Return whether or
        //  not the weight field was found.
        template<template<class> class GeoField>
        bool multiplyWeight
        (
            const word& weightFieldName,
            scalarField& weight
        ) const;

        //- Write the field name for the given field. Return whether or not the
        //  field was found.
        template<template<class> class GeoField, class Type>
        bool writeFieldName(const word& fieldName);

        //- Write the column values for the given field name. Return whether or
        //  not the field was found.
        template<template<class> class GeoField, class Type>
        bool writeFieldValue
        (
            const scalarField& weight,
            const word& fieldName
        );

        //- Write file header information
        virtual void writeFileHeader(const label i = 0);


public:

    //- Runtime type information
    TypeName("LagrangianFieldValue");


    // Constructors

        //- Construct from Time and dictionary
        LagrangianFieldValue
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        LagrangianFieldValue(const LagrangianFieldValue&) = delete;


    //- Destructor
    virtual ~LagrangianFieldValue();


    // Member Functions

        //- Read parameters
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute. Does nothing.
        virtual bool execute();

        //- Write the sum
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LagrangianFieldValue&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
