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
    Foam::interfacialGrowthSizeGroupFvScalarFieldSource

Description
    This source condition creates a interfacial mass transfer source term in
    the size group equations of a population balance.

Usage
    Example of the source condition specification:
    \verbatim
    <sourceName>
    {
        type            interfacialGrowthSizeGroup;
        libs            ("libmultiphaseEulerFvModels.so");
    }
    \endverbatim

SourceFiles
    interfacialGrowthSizeGroupFvScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef interfacialGrowthSizeGroupFvScalarFieldSource_H
#define interfacialGrowthSizeGroupFvScalarFieldSource_H

#include "growthSizeGroupFvScalarFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class interfacialGrowthSizeGroupFvScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class interfacialGrowthSizeGroupFvScalarFieldSource
:
    public growthSizeGroupFvScalarFieldSource
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
    TypeName("interfacialGrowthSizeGroup");


    // Constructors

        //- Inherit constructors
        using
            growthSizeGroupFvScalarFieldSource::
            growthSizeGroupFvScalarFieldSource;

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<fvScalarFieldSource> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return autoPtr<fvScalarFieldSource>
            (
                new interfacialGrowthSizeGroupFvScalarFieldSource
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
