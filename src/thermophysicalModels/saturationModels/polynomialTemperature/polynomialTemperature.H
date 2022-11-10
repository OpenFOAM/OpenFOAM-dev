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
    Foam::saturationModels::polynomialTemperature

Description
    Polynomial equation for the saturation vapour temperature in terms of
    the vapour pressure (in Pa).

    \f[
        T_sat = \sum_i C_i p^i
    \f]

    where \f$p\f$ is the pressure in Pa and \f$C\f$ are the coefficients.

SourceFiles
    polynomialTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef polynomialTemperature_H
#define polynomialTemperature_H

#include "saturationTemperatureModel.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{

/*---------------------------------------------------------------------------*\
                         Class polynomialTemperature Declaration
\*---------------------------------------------------------------------------*/

class polynomialTemperature
:
    public saturationTemperatureModel
{
    // Private Data

        //- Polynomial coefficients
        Polynomial<8> C_;


public:

    //- Runtime type information
    TypeName("polynomial");


    // Constructors

        //- Construct from a dictionary
        polynomialTemperature(const dictionary& dict);


    //- Destructor
    virtual ~polynomialTemperature();


    // Member Functions

        //- Saturation temperature for volScalarField::Internal
        virtual tmp<volScalarField::Internal> Tsat
        (
            const volScalarField::Internal& p
        ) const;

        //- Saturation temperature for volScalarField
        virtual tmp<volScalarField> Tsat(const volScalarField& p) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace saturationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
