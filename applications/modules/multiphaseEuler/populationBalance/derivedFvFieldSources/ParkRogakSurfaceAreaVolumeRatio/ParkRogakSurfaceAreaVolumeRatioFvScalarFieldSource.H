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
    Foam::ParkRogakSurfaceAreaVolumeRatioFvScalarFieldSource

Description
    This source condition creates an interfacial mass transfer source term in
    the surface-area-to-volume-ratio (kappa) equations of the fractal shape
    model of population balance, according to the model of Park and Rogak.

    Reference:
    \verbatim
        Park, S. H., & Rogak, S. N. (2004).
        A novel fixed-sectional model for the formation and growth of aerosol
        agglomerates.
        Journal of Aerosol Science, 35(11), 1385-1404.
    \endbverbatim

Usage
    Example of the source condition specification:
    \verbatim
    <sourceName>
    {
        type            ParkRogakSurfaceAreaVolumeRatio;
        libs            ("libmultiphaseEulerFvModels.so");
    }
    \endverbatim

SourceFiles
    ParkRogakSurfaceAreaVolumeRatioFvScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef ParkRogakSurfaceAreaVolumeRatioFvScalarFieldSource_H
#define ParkRogakSurfaceAreaVolumeRatioFvScalarFieldSource_H

#include "growthSecondaryPropertyFvScalarFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class ParkRogakSurfaceAreaVolumeRatioFvScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class ParkRogakSurfaceAreaVolumeRatioFvScalarFieldSource
:
    public growthSecondaryPropertyFvScalarFieldSource
{
public:

    //- Runtime type information
    TypeName("ParkRogakSurfaceAreaVolumeRatio");


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
                new ParkRogakSurfaceAreaVolumeRatioFvScalarFieldSource
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Return the value transported from the adjacent group
        virtual tmp<DimensionedField<scalar, volMesh>> value
        (
            const label i,
            const fvSource& model
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
