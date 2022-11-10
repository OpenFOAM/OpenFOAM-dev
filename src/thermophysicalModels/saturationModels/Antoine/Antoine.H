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
    Foam::saturationModels::Antoine

Description
    Antoine equation for the vapour pressure.

    \f[
        \log p = A + \frac{B}{C + T}
    \f]

    Coefficients \f$A\f$, \f$B\f$ and \f$C\f$ are to be supplied and should be
    suitable for natural logarithms and temperatures in Kelvin.

SourceFiles
    Antoine.C

\*---------------------------------------------------------------------------*/

#ifndef Antoine_H
#define Antoine_H

#include "saturationPressureModel.H"
#include "saturationTemperatureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{

/*---------------------------------------------------------------------------*\
                           Class Antoine Declaration
\*---------------------------------------------------------------------------*/

class Antoine
:
    public saturationPressureModel,
    public saturationTemperatureModel
{
    // Private data

        //- Constant A
        dimensionedScalar A_;

        //- Constant B
        dimensionedScalar B_;

        //- Constant C
        dimensionedScalar C_;


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

        //- Saturation temperature
        template<class FieldType>
        tmp<FieldType> Tsat(const FieldType& p) const;


public:

    //- Runtime type information
    TypeName("Antoine");


    // Constructors

        //- Construct from a dictionary
        Antoine(const dictionary& dict);


    //- Destructor
    virtual ~Antoine();


    // Member Functions

        //- Saturation pressure for volScalarField::Internal
        DEFINE_PSAT(volScalarField::Internal, );

        //- Saturation pressure for volScalarField
        DEFINE_PSAT(volScalarField, );

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
