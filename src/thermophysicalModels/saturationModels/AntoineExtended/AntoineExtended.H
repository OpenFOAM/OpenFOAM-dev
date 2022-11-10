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
    Foam::saturationModels::AntoineExtended

Description
    Extended Antoine equation for the vapour pressure.

    \f[
        \log (p) = A + \frac{B}{C + T} + D \log (T) + E T^F
    \f]

    Coefficients \f$A\f$, \f$B\f$, \f$C\f$, \f$D\f$, \f$E\f$ and \f$F\f$ are
    to be supplied and should be suitable for natural logarithms and
    temperatures in Kelvin.

SourceFiles
    AntoineExtended.C

\*---------------------------------------------------------------------------*/

#ifndef AntoineExtended_H
#define AntoineExtended_H

#include "saturationPressureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{

/*---------------------------------------------------------------------------*\
                       Class AntoineExtended Declaration
\*---------------------------------------------------------------------------*/

class AntoineExtended
:
    public saturationPressureModel
{
    // Private data

        //- Constant A
        dimensionedScalar A_;

        //- Constant B
        dimensionedScalar B_;

        //- Constant C
        dimensionedScalar C_;

        //- Constant D
        dimensionedScalar D_;

        //- Constant F
        dimensionedScalar F_;

        //- Constant E
        //  (after F so F's dimensions can be used in the construction)
        dimensionedScalar E_;


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
    TypeName("AntoineExtended");


    // Constructors

        //- Construct from a dictionary and an interface
        AntoineExtended(const dictionary& dict);


    //- Destructor
    virtual ~AntoineExtended();


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
