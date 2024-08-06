/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

#include "isotropic.H"
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SolidThermophysicalTransportModel>
Foam::solidThermophysicalTransportModels::
isotropic<SolidThermophysicalTransportModel>::isotropic
(
    const alphaField& alpha,
    const solidThermo& thermo
)
:
    SolidThermophysicalTransportModel(typeName, alpha, thermo)
{
    if (!thermo.isotropic())
    {
        FatalIOErrorInFunction(*this)
            << "Cannot instantiate an isotropic transport model "
               "with anisotropic solid properties"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SolidThermophysicalTransportModel>
bool Foam::solidThermophysicalTransportModels::
isotropic<SolidThermophysicalTransportModel>::read()
{
    return true;
}


template<class SolidThermophysicalTransportModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::solidThermophysicalTransportModels::
isotropic<SolidThermophysicalTransportModel>::q() const
{
    return surfaceScalarField::New
    (
        "q",
       -fvc::interpolate(this->alpha()*this->kappa())
       *fvc::snGrad(this->thermo().T())
    );
}


template<class SolidThermophysicalTransportModel>
Foam::tmp<Foam::scalarField>
Foam::solidThermophysicalTransportModels::
isotropic<SolidThermophysicalTransportModel>::qCorr
(
    const label patchi
) const
{
    return tmp<scalarField>(nullptr);
}


template<class SolidThermophysicalTransportModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::solidThermophysicalTransportModels::
isotropic<SolidThermophysicalTransportModel>::divq
(
    volScalarField& e
) const
{
    const solidThermo& thermo = this->thermo();

    // Return heat flux source as an implicit energy correction
    // to the temperature gradient flux
    return
       -fvc::laplacian(this->alpha()*this->kappa(), thermo.T())
       -fvm::laplacianCorrection(this->alpha()*this->kappa()/thermo.Cv(), e);
}


template<class SolidThermophysicalTransportModel>
void Foam::solidThermophysicalTransportModels::
isotropic<SolidThermophysicalTransportModel>::predict()
{
    SolidThermophysicalTransportModel::predict();
}


// ************************************************************************* //
