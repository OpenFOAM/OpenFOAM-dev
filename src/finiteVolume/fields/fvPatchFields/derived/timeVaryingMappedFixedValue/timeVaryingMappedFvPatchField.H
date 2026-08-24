/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::timeVaryingMappedFvPatchField

Description
    Patch field mapper which interpolates the values from a set of supplied
    points in space and time.

    By default the data files should be provide in
    constant/boundaryData/\<patch name\>/ directory:
      - points             : pointField of locations
      - \<time\>/\<field\> : field of values at time \<time\>

    Alternatively the names and locations of the points and field files may be
    specified explicitly via the optional dictionary entries:
      - dataDir \<optional top-level directory of the points and field data>;
      - points \<optional path including name of points file relative to
                 dataDir\>;
      - sample \<optional name of the sub-directory in the time directories
                containing the fields\>;
    This is particularly useful when mapping data from another case for which
    the \c sample \c functionObject is used to obtain the patch field data for
    mapping.

    If the surface from which the data was obtained changes in time, e.g. when
    sampled from a case with a changing mesh, a separate points file may be
    provided in each time directory alongside the field data, i.e.
    dataDir/\<time\>/\<sample\>/\<points name\>. Points files found in the
    time directories take precedence over the static points file.

    The default mode of operation (mapMethod planarInterpolation) is to project
    the points onto a plane (constructed from the first three points) and
    construct a 2D triangulation and finds for the face centres the triangle it
    is in and the weights to the 3 vertices.

    The optional mapMethod nearest will avoid all projection and triangulation
    and just use the value at the nearest vertex.

    Values are interpolated linearly between times.

Usage
    \table
        Property     | Description                | Required    | Default value
        setAverage   | Switch to activate setting of average value | no | false
        perturb      | Perturb points for regular geometries | no | 1e-5
        fieldTable   | Alternative field name to sample | no| this field name
        mapMethod    | Type of mapping            | no | planarInterpolation
        offset       | Offset to mapped values    | no | Zero
        dataDir      | Top-level directory of the points and field data \\
                         | no | constant/boundaryData/\<patch name\>
        points       | Path including name of points file relative to dataDir \\
                         | no | points
        sample       | Name of the sub-directory in the time directories \\
                         containing the fields | no | ""
    \endtable

See also
    Foam::timeVaryingMappedFixedValueFvPatchField
    Foam::Function1s

SourceFiles
    timeVaryingMappedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingMappedFvPatchField_H
#define timeVaryingMappedFvPatchField_H

#include "fvPatchField.H"
#include "timeVaryingMappedPoints.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

TemplateName(timeVaryingMappedFvPatchField);

/*---------------------------------------------------------------------------*\
                Class timeVaryingMappedFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class timeVaryingMappedFvPatchField
:
    public timeVaryingMappedFvPatchFieldName
{
    // Private Data

        //- Reference to the patch the data is mapped to
        const fvPatch& patch_;

        //- Reference to internal field
        const DimensionedField<Type, fvMesh>& internalField_;

        //- Name of the field data table, defaults to the name of the field
        word fieldTableName_;

        //- If true adjust the mapped field to maintain average value
        Switch setAverage_;

        //- Sample points, times, field files and mapper construction
        timeVaryingMappedPoints points_;

        //- Current starting index in the sample times
        label startSampleTime_;

        //- Interpolated values from startSampleTime
        Field<Type> startSampledValues_;

        //- If setAverage: starting average value
        Type startAverage_;

        //- Current end index in the sample times
        label endSampleTime_;

        //- Interpolated values from endSampleTime
        Field<Type> endSampledValues_;

        //- If setAverage: end average value
        Type endAverage_;

        //- Time varying offset values to interpolated data
        autoPtr<Function1<Type>> offset_;

        //- Read and map the field values for the given sample time index,
        //  setting the average if setAverage is true
        tmp<Field<Type>> sampledValues(const label sampleIndex, Type& average);

        //- Find boundary data in between current time and interpolate
        void checkTable();


public:

    // Constructors

        //- Construct from patch, internal field and dictionary
        timeVaryingMappedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, fvMesh>&,
            const dictionary&
        );

        //- Construct by mapping given timeVaryingMappedFvPatchField
        //  onto a new patch
        timeVaryingMappedFvPatchField
        (
            const timeVaryingMappedFvPatchField&,
            const fvPatch&,
            const DimensionedField<Type, fvMesh>&,
            const fieldMapper&
        );

        //- Copy constructor
        timeVaryingMappedFvPatchField
        (
            const timeVaryingMappedFvPatchField<Type>&
        );


    // Member Functions

        //- Return reference to time
        const Time& time() const
        {
            return patch_.time();
        }

        // Mapping functions

            //- Map the given timeVaryingMappedFvPatchField
            //  onto this timeVaryingMappedFvPatchField
            void map
            (
                const timeVaryingMappedFvPatchField<Type>&,
                const fieldMapper&
            );

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            void reset(const timeVaryingMappedFvPatchField&);


        // Evaluation functions

            //- Return the current mapped patch field
            tmp<Field<Type>> map();


        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "timeVaryingMappedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
