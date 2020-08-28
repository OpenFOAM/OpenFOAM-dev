/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
#include "FixedList.H"
#include "instantList.H"
#include "pointToPointPlanarInterpolation.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class timeVaryingMapped Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingMapped
{
public:

    //- Runtime type information
    ClassName("timeVaryingMapped");


    // Constructors

        //- Construct null
        timeVaryingMapped()
        {}
};


/*---------------------------------------------------------------------------*\
                Class timeVaryingMappedFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class timeVaryingMappedFvPatchField
:
    public timeVaryingMapped
{
    // Private Data

        //- Reference to the patch the data is mapped to
        const fvPatch& patch_;

        //- Name of the field data table, defaults to the name of the field
        word fieldTableName_;

        //- Optional name of the data directory,
        //  defaults to "constant/boundaryData/<patch name>"
        fileName dataDir_;

        //- Optional name of the points file, defaults to "points"
        fileName pointsName_;

        //- Optional name of the sample sub-directory from which the field are
        //  read, defaults to ""
        fileName sampleName_;

        //- If true adjust the mapped field to maintain average value
        Switch setAverage_;

        //- Fraction of perturbation (fraction of bounding box) to add
        scalar perturb_;

        //- Interpolation scheme to use
        word mapMethod_;

        //- 2D interpolation (for 'planarInterpolation' mapMethod)
        autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

        //- List of boundaryData time directories
        instantList sampleTimes_;

        //- Current starting index in sampleTimes
        label startSampleTime_;

        //- Interpolated values from startSampleTime
        Field<Type> startSampledValues_;

        //- If setAverage: starting average value
        Type startAverage_;

        //- Current end index in sampleTimes
        label endSampleTime_;

        //- Interpolated values from endSampleTime
        Field<Type> endSampledValues_;

        //- If setAverage: end average value
        Type endAverage_;

        //- Time varying offset values to interpolated data
        autoPtr<Function1<Type>> offset_;

        //- Helper function to find the field files
        fileName findFieldFile(const word& timeName) const;

        //- Find boundary data in between current time and interpolate
        void checkTable();


public:

    // Constructors

        //- Construct from patch and field name
        timeVaryingMappedFvPatchField(const fvPatch&, const word& fieldName);

        //- Construct from patch, dictionary and field name
        timeVaryingMappedFvPatchField
        (
            const fvPatch&,
            const dictionary&,
            const word& fieldName
        );

        //- Copy constructor
        timeVaryingMappedFvPatchField
        (
            const timeVaryingMappedFvPatchField<Type>&
        );


    // Member Functions

        const Time& time() const
        {
            return patch_.boundaryMesh().mesh().time();
        }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given timeVaryingMappedFvPatchField
            //  onto this timeVaryingMappedFvPatchField
            void rmap
            (
                const timeVaryingMappedFvPatchField<Type>&,
                const labelList&
            );


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
