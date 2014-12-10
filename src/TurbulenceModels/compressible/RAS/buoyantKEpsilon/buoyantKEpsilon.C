/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "buoyantKEpsilon.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
buoyantKEpsilon<BasicTurbulenceModel>::buoyantKEpsilon
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kEpsilon<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    Cg_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cg",
            this->coeffDict_,
            0.85
        )
    )
{
    if (type == typeName)
    {
        kEpsilon<BasicTurbulenceModel>::correctNut();
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool buoyantKEpsilon<BasicTurbulenceModel>::read()
{
    if (kEpsilon<BasicTurbulenceModel>::read())
    {
        Cg_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField>
buoyantKEpsilon<BasicTurbulenceModel>::Gcoef() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    return
        this->alpha_*(this->Cmu_/this->sigmak_)*this->k_
       *(g & fvc::grad(this->rho_))/(this->epsilon_ + this->epsilonMin_);
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
buoyantKEpsilon<BasicTurbulenceModel>::kSource() const
{
    return -fvm::SuSp(Gcoef(), this->k_);
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
buoyantKEpsilon<BasicTurbulenceModel>::epsilonSource() const
{
    return -fvm::SuSp(this->C1_*(1.0 - this->Cg_)*Gcoef(), this->epsilon_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
