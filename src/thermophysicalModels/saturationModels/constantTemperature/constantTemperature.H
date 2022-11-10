/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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
    Foam::saturationModels::constantTemperature

Description
    Constant saturation temperature model.

SourceFiles
    constantTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef constantTemperature_H
#define constantTemperature_H

#include "saturationTemperatureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{

/*---------------------------------------------------------------------------*\
                Class constantTemperature Declaration
\*---------------------------------------------------------------------------*/

class constantTemperature
:
    public saturationTemperatureModel
{
protected:

    // Private Data

        //- Constant saturation temperature
        const dimensionedScalar Tsat_;


    // Private Member Functions

        //- Saturation temperature
        template<class FieldType>
        tmp<FieldType> Tsat(const FieldType& p) const;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from a dictionary
        constantTemperature(const dictionary& dict);

        //- Construct from a value
        constantTemperature(const dimensionedScalar& Tsat);


    //- Destructor
    virtual ~constantTemperature();


    // Member Functions

        //- Saturation pressure for volScalarField::Internal
        DEFINE_TSAT(volScalarField::Internal, );

        //- Saturation pressure for volScalarField
        DEFINE_TSAT(volScalarField, );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace saturationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
