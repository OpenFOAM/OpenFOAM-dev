/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::energyCalculatedTemperatureFvScalarFieldSource

Description
    Base class for temperature source conditions in which the parameters of the
    corresponding energy condition can be set directly.

SourceFiles
    energyCalculatedTemperatureFvScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef energyCalculatedTemperatureFvScalarFieldSource_H
#define energyCalculatedTemperatureFvScalarFieldSource_H

#include "fvFieldSources.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class energyCalculatedTemperatureFvScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class energyCalculatedTemperatureFvScalarFieldSource
:
    public fvScalarFieldSource
{
public:

    //- Runtime type information
    TypeName("energyCalculatedTemperature");


    // Constructors

        //- Construct from internal field and dictionary
        energyCalculatedTemperatureFvScalarFieldSource
        (
            const DimensionedField<scalar, volMesh>&,
            const dictionary& dict
        );

        //- Copy constructor setting internal field reference
        energyCalculatedTemperatureFvScalarFieldSource
        (
            const energyCalculatedTemperatureFvScalarFieldSource&,
            const DimensionedField<scalar, volMesh>&
        );


    //- Destructor
    virtual ~energyCalculatedTemperatureFvScalarFieldSource();


    // Member Functions

        //- Return the source energy value
        virtual tmp<scalarField> sourceHeValue(const fvSource&) const = 0;

        //- Return the source value. Not implemented.
        virtual tmp<scalarField> sourceValue(const fvSource&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
