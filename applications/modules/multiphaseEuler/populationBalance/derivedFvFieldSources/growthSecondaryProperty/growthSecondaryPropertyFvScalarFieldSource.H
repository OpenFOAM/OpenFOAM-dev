/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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
    Foam::growthSecondaryPropertyFvScalarFieldSource

Description
    Base class for source conditions which create a growth source terms for
    secondary properties of the size-groups of a population balance model.

SourceFiles
    growthSecondaryPropertyFvScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef growthSecondaryPropertyFvScalarFieldSource_H
#define growthSecondaryPropertyFvScalarFieldSource_H

#include "growthFvScalarFieldSource.H"
#include "secondaryPropertyFvScalarFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class growthSecondaryPropertyFvScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class growthSecondaryPropertyFvScalarFieldSource
:
    public growthFvScalarFieldSource,
    public secondaryPropertyFvScalarFieldSource
{
public:

    // Constructors

        //- Construct from internal field and dictionary
        growthSecondaryPropertyFvScalarFieldSource
        (
            const DimensionedField<scalar, volMesh>&,
            const dictionary& dict
        );

        //- Copy constructor setting internal field reference
        growthSecondaryPropertyFvScalarFieldSource
        (
            const growthSecondaryPropertyFvScalarFieldSource&,
            const DimensionedField<scalar, volMesh>&
        );


    // Member Functions

        //- Return the value transported from the adjacent group
        virtual tmp<DimensionedField<scalar, volMesh>> value
        (
            const label i,
            const fvSource& model
        ) const = 0;

        //- Return the internal coefficient
        virtual tmp<DimensionedField<scalar, volMesh>> internalCoeff
        (
            const fvSource& model,
            const DimensionedField<scalar, volMesh>& source
        ) const;

        //- Return the source coefficients for exchange with the groups below
        //  and above this group
        virtual Pair<tmp<DimensionedField<scalar, volMesh>>> sourceCoeffs
        (
            const fvSource& model
        ) const;

        //- Return the combined source coefficient
        virtual tmp<DimensionedField<scalar, volMesh>> sourceCoeff
        (
            const fvSource& model,
            const DimensionedField<scalar, volMesh>& source
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
