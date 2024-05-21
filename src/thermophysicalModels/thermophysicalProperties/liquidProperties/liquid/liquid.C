/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2024 OpenFOAM Foundation
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

#include "liquid.H"
#include "None.H"
#include "addToRunTimeSelectionTable.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(liquid, 0);
    addToRunTimeSelectionTable(liquidProperties, liquid, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::Function1<Foam::scalar>> Foam::liquid::New
(
    const word& name,
    const dimensionSet& dims,
    const dictionary& dict
)
{
    if (dict.isDict(name))
    {
        return Function1<scalar>::New(name, dimTemperature, dims, dict);
    }
    else
    {
        return autoPtr<Function1<scalar>>(new Function1s::None<scalar>(name));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquid::liquid(const dictionary& dict)
:
    liquidProperties(dict),
    rho_(New("rho", dimDensity, dict)),
    pv_(New("pv", dimPressure, dict)),
    hl_(New("hl", dimEnergy/dimMass, dict)),
    Cp_(New("Cp", dimSpecificHeatCapacity, dict)),
    h_(New("h", dimEnergy/dimMass, dict)),
    Cpg_(New("Cpg", dimSpecificHeatCapacity, dict)),
    B_(New("B", dimVolume/dimMass, dict)),
    mu_(New("mu", dimDynamicViscosity, dict)),
    mug_(New("mug", dimDynamicViscosity, dict)),
    kappa_(New("kappa", dimThermalConductivity, dict)),
    kappag_(New("kappag", dimThermalConductivity, dict)),
    sigma_(New("sigma", dimForce/dimLength, dict)),
    D_(New("D", dimArea/dimTime, dict)),
    hf_(h_->value(Tstd))
{}


Foam::liquid::liquid(const liquid& lm)
:
    liquidProperties(lm),
    rho_(lm.rho_, false),
    pv_(lm.pv_, false),
    hl_(lm.hl_, false),
    Cp_(lm.Cp_, false),
    h_(lm.h_, false),
    Cpg_(lm.Cpg_, false),
    B_(lm.B_, false),
    mu_(lm.mu_, false),
    mug_(lm.mug_, false),
    kappa_(lm.kappa_, false),
    kappag_(lm.kappag_, false),
    sigma_(lm.sigma_, false),
    D_(lm.D_, false),
    hf_(lm.hf_)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::liquid::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
