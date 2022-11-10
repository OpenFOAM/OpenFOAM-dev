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
    Foam::saturationModels::constantPressure

Description
    Constant saturation pressure model.

SourceFiles
    constantPressure.C

\*---------------------------------------------------------------------------*/

#ifndef constantPressure_H
#define constantPressure_H

#include "saturationPressureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{

/*---------------------------------------------------------------------------*\
                          Class constantPressure Declaration
\*---------------------------------------------------------------------------*/

class constantPressure
:
    public saturationPressureModel
{
    // Private Data

        //- The constant saturation pressure
        const dimensionedScalar pSat_;


    // Private Member Functions

        //- Saturation pressure
        template<class FieldType>
        tmp<FieldType> pSat(const FieldType& T) const;

        //- Saturation pressure derivative w.r.t. temperature
        template<class FieldType>
        tmp<FieldType> pSatPrime(const FieldType& T) const;

        //- Natural log of the saturation pressure
        template<class FieldType>
        tmp<FieldType> lnPSat(const FieldType& T) const;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from a dictionary
        constantPressure(const dictionary& dict);

        //- Construct from a value
        constantPressure(const dimensionedScalar& pSat);


    //- Destructor
    virtual ~constantPressure();


    // Member Functions

        //- Saturation pressure for volScalarField::Internal
        DEFINE_PSAT(volScalarField::Internal, );

        //- Saturation pressure for volScalarField
        DEFINE_PSAT(volScalarField, );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace saturationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
