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
    Foam::interfacialGrowthGroupFractionFvScalarFieldSource

Description
    This source condition creates an interfacial growth transfer source term in
    the group fraction equations of a population balance model, consistent with
    interfacial mass transfer of phase change or similar.

Usage
    Example of the source condition specification:
    \verbatim
    <sourceName>
    {
        type            interfacialGrowthGroupFraction;
        libs            ("libmultiphaseEulerFvModels.so");
    }
    \endverbatim

SourceFiles
    interfacialGrowthGroupFractionFvScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef interfacialGrowthGroupFractionFvScalarFieldSource_H
#define interfacialGrowthGroupFractionFvScalarFieldSource_H

#include "growthGroupFractionFvScalarFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class interfacialGrowthGroupFractionFvScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class interfacialGrowthGroupFractionFvScalarFieldSource
:
    public growthGroupFractionFvScalarFieldSource
{
private:

    // Private Member Functions

        //- Return the moment of the diameter with which the source term scales
        virtual label q() const
        {
            return 2;
        }


public:

    //- Runtime type information
    TypeName("interfacialGrowthGroupFraction");


    // Constructors

        //- Inherit constructors
        using
            growthGroupFractionFvScalarFieldSource::
            growthGroupFractionFvScalarFieldSource;

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<fvScalarFieldSource> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return autoPtr<fvScalarFieldSource>
            (
                new interfacialGrowthGroupFractionFvScalarFieldSource
                (
                    *this,
                    iF
                )
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
