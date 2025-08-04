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
    Foam::interfacialGrowthGroupSurfaceAreaVolumeRatioFvScalarFieldSource

Description
    This source condition creates an interfacial growth transfer source term in
    the group surface-area-to-volume-ratio (kappa) equations of the fractal
    shape model of a population balance model, consistent with interfacial mass
    transfer of phase change or similar.

Usage
    Example of the source condition specification:
    \verbatim
    <sourceName>
    {
        type            interfacialGrowthGroupSurfaceAreaVolumeRatio;
        libs            ("libmultiphaseEulerFvModels.so");
    }
    \endverbatim

SourceFiles
    interfacialGrowthGroupSurfaceAreaVolumeRatioFvScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef interfacialGrowthGroupSurfaceAreaVolumeRatioFvScalarFieldSource_H
#define interfacialGrowthGroupSurfaceAreaVolumeRatioFvScalarFieldSource_H

#include "growthSecondaryPropertyFvScalarFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
Class interfacialGrowthGroupSurfaceAreaVolumeRatioFvScalarFieldS... Declaration
\*---------------------------------------------------------------------------*/

class interfacialGrowthGroupSurfaceAreaVolumeRatioFvScalarFieldSource
:
    public growthSecondaryPropertyFvScalarFieldSource
{
public:

    //- Runtime type information
    TypeName("interfacialGrowthGroupSurfaceAreaVolumeRatio");


    // Constructors

        //- Inherit constructors
        using
            growthSecondaryPropertyFvScalarFieldSource::
            growthSecondaryPropertyFvScalarFieldSource;

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<fvScalarFieldSource> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return autoPtr<fvScalarFieldSource>
            (
                new
                interfacialGrowthGroupSurfaceAreaVolumeRatioFvScalarFieldSource
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Return the value transported from an adjacent group
        virtual tmp<DimensionedField<scalar, volMesh>> value
        (
            const label j,
            const fvSource& model
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
