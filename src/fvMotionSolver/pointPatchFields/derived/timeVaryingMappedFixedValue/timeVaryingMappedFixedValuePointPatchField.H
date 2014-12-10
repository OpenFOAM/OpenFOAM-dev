/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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
    Foam::timeVaryingMappedFixedValuePointPatchField

Description
    A time-varying form of a mapped fixed value boundary condition.

See Also
    Foam::timeVaryingMappedFixedValueFvPatchField

SourceFiles
    timeVaryingMappedFixedValuePointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingMappedFixedValuePointPatchField_H
#define timeVaryingMappedFixedValuePointPatchField_H

#include "fixedValuePointPatchField.H"
#include "instantList.H"
#include "pointToPointPlanarInterpolation.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class timeVaryingMappedFixedValuePointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class timeVaryingMappedFixedValuePointPatchField
:
    public fixedValuePointPatchField<Type>
{
    // Private data

        //- Name of the field data table, defaults to the name of the field
        word fieldTableName_;

        //- If true adjust the mapped field to maintain average value
        bool setAverage_;

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
        autoPtr<DataEntry<Type> > offset_;


public:

    //- Runtime type information
    TypeName("timeVaryingMappedFixedValue");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingMappedFixedValuePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingMappedFixedValuePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        timeVaryingMappedFixedValuePointPatchField
        (
            const timeVaryingMappedFixedValuePointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingMappedFixedValuePointPatchField
        (
            const timeVaryingMappedFixedValuePointPatchField<Type>&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<Type> > clone() const
        {
            return autoPtr<pointPatchField<Type> >
            (
                new timeVaryingMappedFixedValuePointPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingMappedFixedValuePointPatchField
        (
            const timeVaryingMappedFixedValuePointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type> > clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<Type> >
            (
                new timeVaryingMappedFixedValuePointPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Utility functions

            //- Find boundary data inbetween current time and interpolate
            void checkTable();


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const pointPatchFieldMapper&
            );

            //- Reverse map the given PointPatchField onto
            //  this PointPatchField
            virtual void rmap
            (
                const pointPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "timeVaryingMappedFixedValuePointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
