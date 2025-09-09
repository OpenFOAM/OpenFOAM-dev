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
    Foam::interfacialGrowthInterfacialCurvatureFvScalarFieldSource

Description
    This source condition creates a interfacial mass transfer source term in
    the interfacial curvature equations of the IATE model.

Usage
    Example of the source condition specification:
    \verbatim
    <sourceName>
    {
        type            interfacialGrowthInterfacialCurvature;
        libs            ("libmultiphaseEulerFvModels.so");
    }
    \endverbatim

SourceFiles
    interfacialGrowthInterfacialCurvatureFvScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef interfacialGrowthInterfacialCurvatureFvScalarFieldSource_H
#define interfacialGrowthInterfacialCurvatureFvScalarFieldSource_H

#include "fvFieldSources.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class interfacialGrowthInterfacialCurvatureFvScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class interfacialGrowthInterfacialCurvatureFvScalarFieldSource
:
    public fvScalarFieldSource
{
public:

    //- Runtime type information
    TypeName("interfacialGrowthInterfacialCurvature");


    // Constructors

        //- Construct from internal field and dictionary
        interfacialGrowthInterfacialCurvatureFvScalarFieldSource
        (
            const DimensionedField<scalar, volMesh>&,
            const dictionary& dict
        );

        //- Copy constructor setting internal field reference
        interfacialGrowthInterfacialCurvatureFvScalarFieldSource
        (
            const interfacialGrowthInterfacialCurvatureFvScalarFieldSource&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<fvScalarFieldSource> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return autoPtr<fvScalarFieldSource>
            (
                new interfacialGrowthInterfacialCurvatureFvScalarFieldSource
                (
                    *this,
                    iF
                )
            );
        }


    //- Destructor
    virtual ~interfacialGrowthInterfacialCurvatureFvScalarFieldSource();


    // Member Functions

        //- Return the source value
        virtual tmp<DimensionedField<scalar, volMesh>> sourceValue
        (
            const fvSource& model,
            const DimensionedField<scalar, volMesh>& source
        ) const;

        //- Return the internal coefficient
        virtual tmp<DimensionedField<scalar, volMesh>> internalCoeff
        (
            const fvSource& model,
            const DimensionedField<scalar, volMesh>& source
        ) const;

        //- Return the source coefficient
        virtual tmp<DimensionedField<scalar, volMesh>> sourceCoeff
        (
            const fvSource& model,
            const DimensionedField<scalar, volMesh>& source
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
