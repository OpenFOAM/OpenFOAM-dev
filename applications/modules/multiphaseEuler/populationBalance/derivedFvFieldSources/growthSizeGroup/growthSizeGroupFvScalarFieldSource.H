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
    Foam::growthSizeGroupFvScalarFieldSource

Description
    Base class for source conditions which create a growth source terms in the
    size group equations of a population balance.

SourceFiles
    growthSizeGroupFvScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef growthSizeGroupFvScalarFieldSource_H
#define growthSizeGroupFvScalarFieldSource_H

#include "growthFvScalarFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
namespace diameterModels
{
    class sizeGroup;
}

/*---------------------------------------------------------------------------*\
              Class growthSizeGroupFvScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class growthSizeGroupFvScalarFieldSource
:
    public growthFvScalarFieldSource
{
private:

    // Private Member Functions

        //- Return the moment of the diameter with which the source term scales
        virtual label q() const = 0;

        //- Compute and return the weight for the given size-group
        tmp<DimensionedField<scalar, volMesh>> w
        (
            const fvSource&,
            const diameterModels::sizeGroup&
        ) const;


public:

    // Constructors

        //- Inherit constructors
        using growthFvScalarFieldSource::growthFvScalarFieldSource;


    // Member Functions

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
