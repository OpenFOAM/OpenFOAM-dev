/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::angularOscillatingVelocityPointPatchVectorField

Description
    Foam::angularOscillatingVelocityPointPatchVectorField

SourceFiles
    angularOscillatingVelocityPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef angularOscillatingVelocityPointPatchVectorField_H
#define angularOscillatingVelocityPointPatchVectorField_H

#include "pointPatchFields.H"
#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class angularOscillatingVelocityPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class angularOscillatingVelocityPointPatchVectorField
:
    public fixedValuePointPatchVectorField
{
    // Private Data

        vector axis_;
        vector origin_;
        scalar angle0_;
        scalar amplitude_;
        scalar omega_;

        pointField p0_;


public:

    //- Runtime type information
    TypeName("angularOscillatingVelocity");


    // Constructors

        //- Construct from patch, internal field and dictionary
        angularOscillatingVelocityPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  angularOscillatingVelocityPointPatchVectorField
        //  onto a new patch
        angularOscillatingVelocityPointPatchVectorField
        (
            const angularOscillatingVelocityPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const fieldMapper&
        );

        //- Construct as copy setting internal field reference
        angularOscillatingVelocityPointPatchVectorField
        (
            const angularOscillatingVelocityPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchVectorField> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchVectorField>
            (
                new angularOscillatingVelocityPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Mapping functions

            //- Map the given pointPatchField onto this pointPatchField
            virtual void map(const pointPatchVectorField&, const fieldMapper&);

            //- Reset the pointPatchField to the given pointPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const pointPatchVectorField&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using fixedValuePointPatchVectorField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
