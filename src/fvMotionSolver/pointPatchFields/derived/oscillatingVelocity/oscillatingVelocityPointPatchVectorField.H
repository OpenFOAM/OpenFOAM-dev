/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::oscillatingVelocityPointPatchVectorField

Description
    Foam::oscillatingVelocityPointPatchVectorField

SourceFiles
    oscillatingVelocityPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef oscillatingVelocityPointPatchVectorField_H
#define oscillatingVelocityPointPatchVectorField_H

#include "pointPatchFields.H"
#include "fixedValuePointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class oscillatingVelocityPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class oscillatingVelocityPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private Data

        vector amplitude_;
        scalar omega_;

        pointField p0_;


public:

    //- Runtime type information
    TypeName("oscillatingVelocity");


    // Constructors

        //- Construct from patch, internal field and dictionary
        oscillatingVelocityPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  oscillatingVelocityPointPatchVectorField
        //  onto a new patch
        oscillatingVelocityPointPatchVectorField
        (
            const oscillatingVelocityPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const fieldMapper&
        );

        //- Construct as copy setting internal field reference
        oscillatingVelocityPointPatchVectorField
        (
            const oscillatingVelocityPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new oscillatingVelocityPointPatchVectorField
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
            virtual void reset(const pointPatchField<vector>&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
