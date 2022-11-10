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
    Foam::saturationModels::ArdenBuck

Description
    ArdenBuck equation for the vapour pressure of moist air.

SourceFiles
    ArdenBuck.C

\*---------------------------------------------------------------------------*/

#ifndef ArdenBuck_H
#define ArdenBuck_H

#include "saturationPressureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{

/*---------------------------------------------------------------------------*\
                          Class ArdenBuck Declaration
\*---------------------------------------------------------------------------*/

class ArdenBuck
:
    public saturationPressureModel
{
    // Private Member Functions

        //- Exponent divided by the temperature
        template<class FieldType>
        tmp<FieldType> xByTC(const FieldType& TC) const;

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
    TypeName("ArdenBuck");


    // Constructors

        //- Construct from a dictionary
        ArdenBuck(const dictionary& dict);


    //- Destructor
    virtual ~ArdenBuck();


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
