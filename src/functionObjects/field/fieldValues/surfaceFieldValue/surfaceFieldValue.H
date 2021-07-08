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
    Foam::functionObjects::fieldValues::surfaceFieldValue

Description
    Provides a 'face regionType' variant of the fieldValues function object.

    Given a list of user-specified fields and a selection of mesh (or general
    surface) faces, a number of operations can be performed, such as sums,
    averages and integrations.

    For example, to calculate the volumetric or mass flux across a patch,
    apply the 'orientedSum' operator to the flux field (typically \c phi)

    Examples of function object specification:
    \verbatim
    movingWallPatch
    {
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");

        log             true;
        writeControl    writeTime;
        writeFields     true;

        regionType      patch;
        name            movingWall;

        operation       areaAverage;

        fields
        (
            p
            phi
            U
        );
    }

    surfaceFieldValue1
    {
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");

        log             true;
        writeControl    writeTime;
        writeFields     true;

        surfaceFormat   none;
        regionType      faceZone;
        name            f0;

        operation       sum;

        weightField     alpha1;

        fields
        (
            p
            phi
            U
        );
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: surfaceFieldValue   | yes         |
        log          | write data to standard output | no    | no
        writeFields  | Write the region field values  | yes     |
        writeArea    | Write the area of the surfaceFieldValue | no |
        surfaceFormat | output value format    | no          |
        regionType   | face regionType: see below  | yes         |
        name         | name of face regionType if required  | no |
        operation    | operation to perform    | yes         |
        weightField  | name of field to apply weighting | no |
        weightFields | Names of fields to apply weighting | no |
        scaleFactor  | scale factor            | no          | 1
        fields       | list of fields to operate on | yes    |
    \endtable

    Where \c regionType is defined by
    \plaintable
        faceZone     | requires a 'name' entry to specify the faceZone
        patch        | requires a 'name' entry to specify the patch
        sampledSurface | requires a 'sampledSurfaceDict' sub-dictionary
    \endplaintable

    The \c operation is one of:
    \plaintable
       none          | no operation
       sum           | sum
       sumMag        | sum of component magnitudes
       sumDirection  | sum values which are positive in given direction
       sumDirectionBalance | sum of balance of values in given direction
       orientedSum   | sum with face orientations
       average       | ensemble average
       areaAverage   | area weighted average
       areaIntegrate | area integral
       min           | minimum
       max           | maximum
       minMag        | minimum magnitude
       maxMag        | maximum magnitude
       CoV           | coefficient of variation: standard deviation/mean
       areaNormalAverage | area weighted average in face normal direction
       areaNormalIntegrate | area weighted integral in face normal direction
    \endplaintable

Note
    - Faces on empty patches get ignored.
    - The `oriented' operations will flip the sign of the field so that all the
      normals point in a consistent direction. This is only of relevance when
      summing mesh-oriented fields, such as the flux, on faceZones.
    - If the field is a volField then a \c faceZone can only consist of
      boundary faces, because only these faces have a value associated with
      them. No cell-to-face interpolation is performed.
    - If the field is a surfaceField then the region cannot be a \c
      sampledSurface
    - If a sampledSurface has interpolation set to false then the surface
      face values will be taken directly from the cell that contains the
      surface face centre
    - If a \c sampledSurface has interpolation set to true then the field
      will be interpolated to the vertices, then averaged onto the surface
      faces

See also
    Foam::fieldValues
    Foam::functionObject

SourceFiles
    surfaceFieldValue.C
    surfaceFieldValueTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_surfaceFieldValue_H
#define functionObjects_surfaceFieldValue_H

#include "fieldValue.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class sampledSurface;
class surfaceWriter;

namespace functionObjects
{
namespace fieldValues
{

/*---------------------------------------------------------------------------*\
                      Class surfaceFieldValue Declaration
\*---------------------------------------------------------------------------*/

class surfaceFieldValue
:
    public fieldValue
{

public:

    // Public data types

        //- region type enumeration
        enum class regionTypes
        {
            faceZone,
            patch,
            sampledSurface
        };

        //- region type names
        static const NamedEnum<regionTypes, 3> regionTypeNames_;


        //- Operation type enumeration
        enum class operationType
        {
            none,
            sum,
            sumMag,
            sumDirection,
            sumDirectionBalance,
            orientedSum,
            average,
            areaAverage,
            areaIntegrate,
            min,
            max,
            minMag,
            maxMag,
            CoV,
            areaNormalAverage,
            areaNormalIntegrate
        };

        //- Operation type names
        static const NamedEnum<operationType, 16> operationTypeNames_;


private:

    // Private Member Functions

        //- Set faces to evaluate based on a face zone
        void setFaceZoneFaces();

        //- Set faces to evaluate based on a patch
        void setPatchFaces();

        //- Set faces according to sampledSurface
        void sampledSurfaceFaces(const dictionary&);

        //- Combine mesh faces and points from multiple processors
        void combineMeshGeometry
        (
            faceList& faces,
            pointField& points
        ) const;

        //- Combine surface faces and points from multiple processors
        void combineSurfaceGeometry
        (
            faceList& faces,
            pointField& points
        ) const;

        //- Calculate and return total area of the surfaceFieldValue: sum(magSf)
        scalar totalArea() const;


protected:

    // Protected data

        //- Surface writer
        autoPtr<surfaceWriter> surfaceWriterPtr_;

        //- region type
        regionTypes regionType_;

        //- Operation to apply to values
        operationType operation_;

        //- Weight field names - optional
        wordList weightFieldNames_;

        //- Scale factor - optional
        scalar scaleFactor_;

        //- Total area of the surfaceFieldValue
        scalar totalArea_;

        //- Optionally write the area of the surfaceFieldValue
        bool writeArea_;

        //- Global number of faces
        label nFaces_;


        // If operating on mesh faces (faceZone, patch)

            //- Local list of face IDs
            labelList faceId_;

            //- Local list of patch ID per face
            labelList facePatchId_;

            //- List of +1/-1 representing face flip map
            //  (1 use as is, -1 negate)
            labelList faceSign_;


        // If operating on sampledSurface

            //- Underlying sampledSurface
            autoPtr<sampledSurface> surfacePtr_;


    // Protected Member Functions

        //- Initialise, e.g. face addressing
        void initialise(const dictionary& dict);

        //- Return true if the field name is valid
        template<class Type>
        bool validField(const word& fieldName) const;

        //- Return field values by looking up field name
        template<class Type>
        tmp<Field<Type>> getFieldValues(const word& fieldName) const;

        //- Apply the operation to the values, and return true if successful.
        //  Does nothing unless overloaded below.
        template<class Type, class ResultType>
        bool processValues
        (
            const Field<Type>& values,
            const scalarField& signs,
            const scalarField& weights,
            const vectorField& Sf,
            ResultType& result
        ) const;

        //- Apply Type -> Type operation to the values. Calls
        //  processValuesTypeType.
        template<class Type>
        bool processValues
        (
            const Field<Type>& values,
            const scalarField& signs,
            const scalarField& weights,
            const vectorField& Sf,
            Type& result
        ) const;

        //- Apply Type -> scalar operation to the values
        template<class Type>
        bool processValues
        (
            const Field<Type>& values,
            const scalarField& signs,
            const scalarField& weights,
            const vectorField& Sf,
            scalar& result
        ) const;

        //- Apply scalar -> scalar operation to the values. Tries to apply
        //  scalar -> scalar specific operations, otherwise calls
        //  processValuesTypeType.
        bool processValues
        (
            const Field<scalar>& values,
            const scalarField& signs,
            const scalarField& weights,
            const vectorField& Sf,
            scalar& result
        ) const;

        //- Apply vector -> vector operation to the values.
        bool processValues
        (
            const Field<vector>& values,
            const scalarField& signs,
            const scalarField& weights,
            const vectorField& Sf,
            scalar& result
        ) const;

        //- Apply vector -> vector operation to the values. Tries to apply
        //  vector -> vector specific operations, otherwise calls
        //  processValuesTypeType.
        bool processValues
        (
            const Field<vector>& values,
            const scalarField& signs,
            const scalarField& weights,
            const vectorField& Sf,
            vector& result
        ) const;

        //- Apply a Type -> Type operation to the values
        template<class Type>
        bool processValuesTypeType
        (
            const Field<Type>& values,
            const scalarField& signs,
            const scalarField& weights,
            const vectorField& Sf,
            Type& result
        ) const;

        //- Output file header information
        virtual void writeFileHeader(const label i);


public:

    //- Run-time type information
    TypeName("surfaceFieldValue");


    // Constructors

        //- Construct from name, Time and dictionary
        surfaceFieldValue
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from name, objectRegistry and dictionary
        surfaceFieldValue
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfaceFieldValue();


    // Public Member Functions

        //- Return the region type
        inline const regionTypes& regionType() const;

        //- Return the local list of face IDs
        inline const labelList& faceId() const;

        //- Return the local list of patch ID per face
        inline const labelList& facePatch() const;

        //- Return the list of +1/-1 representing face flip map
        inline const labelList& faceSign() const;

        //- Return the output directory
        inline fileName outputDir() const;

        //- Templated helper function to output field values
        template<class Type>
        bool writeValues
        (
            const word& fieldName,
            const scalarField& signs,
            const scalarField& weights,
            const vectorField& Sf
        );

        //- Templated helper function to output field values
        template<class Type, class ResultType>
        bool writeValues
        (
            const word& fieldName,
            const Field<Type>& values,
            const scalarField& signs,
            const scalarField& weights,
            const vectorField& Sf
        );

        //- Filter a surface field according to faceIds
        template<class Type>
        tmp<Field<Type>> filterField
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& field
        ) const;

        //- Filter a volume field according to faceIds
        template<class Type>
        tmp<Field<Type>> filterField
        (
            const GeometricField<Type, fvPatchField, volMesh>& field
        ) const;

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

#include "surfaceFieldValueI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfaceFieldValueTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
