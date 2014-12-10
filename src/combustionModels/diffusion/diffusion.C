/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "diffusion.H"
#include "fvcGrad.H"

namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
diffusion<CombThermoType, ThermoType>::diffusion
(
    const word& modelType, const fvMesh& mesh
)
:
    singleStepCombustion<CombThermoType, ThermoType>(modelType, mesh),
    C_(readScalar(this->coeffs().lookup("C"))),
    oxidantName_(this->coeffs().template lookupOrDefault<word>("oxidant", "O2"))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
diffusion<CombThermoType, ThermoType>::~diffusion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
void diffusion<CombThermoType, ThermoType>::correct()
{
    this->wFuel_ ==
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0);

    if (this->active())
    {
        this->singleMixturePtr_->fresCorrect();

        const label fuelI = this->singleMixturePtr_->fuelIndex();

        const volScalarField& YFuel =
            this->thermoPtr_->composition().Y()[fuelI];

        if (this->thermoPtr_->composition().contains(oxidantName_))
        {
            const volScalarField& YO2 =
                this->thermoPtr_->composition().Y(oxidantName_);

            this->wFuel_ ==
                C_*this->turbulence().muEff()
               *mag(fvc::grad(YFuel) & fvc::grad(YO2))
               *pos(YFuel)*pos(YO2);
        }
    }
}


template<class CombThermoType, class ThermoType>
bool diffusion<CombThermoType, ThermoType>::read()
{
    if (singleStepCombustion<CombThermoType, ThermoType>::read())
    {
        this->coeffs().lookup("C") >> C_ ;
        this->coeffs().readIfPresent("oxidant", oxidantName_);
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
