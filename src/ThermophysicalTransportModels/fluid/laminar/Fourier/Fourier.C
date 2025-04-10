/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2024 OpenFOAM Foundation
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

#include "Fourier.H"
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
Fourier<BasicThermophysicalTransportModel>::Fourier
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    laminarThermophysicalTransportModel<BasicThermophysicalTransportModel>
    (
        typeName,
        momentumTransport,
        thermo
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
bool Fourier<BasicThermophysicalTransportModel>::read()
{
    return true;
}


template<class TurbulenceThermophysicalTransportModel>
tmp<volScalarField>
Fourier<TurbulenceThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi
) const
{
    FatalErrorInFunction
        << type() << " supports single component systems only, " << nl
        << "    for multi-component transport select"
           " unityLewisFourier"
        << exit(FatalError);

    return tmp<volScalarField>(nullptr);
}


template<class TurbulenceThermophysicalTransportModel>
tmp<scalarField>
Fourier<TurbulenceThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{
    FatalErrorInFunction
        << type() << " supports single component systems only, " << nl
        << "    for multi-component transport select"
           " unityLewisFourier"
        << exit(FatalError);

    return tmp<scalarField>(nullptr);
}


template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField> Fourier<BasicThermophysicalTransportModel>::q() const
{
    const thermoModel& thermo = this->thermo();

    return surfaceScalarField::New
    (
        IOobject::groupName
        (
            "q",
            this->momentumTransport().alphaRhoPhi().group()
        ),
        -fvc::interpolate(this->alpha()*thermo.kappa())*fvc::snGrad(thermo.T())
    );
}


template<class BasicThermophysicalTransportModel>
tmp<scalarField> Fourier<BasicThermophysicalTransportModel>::q
(
    const label patchi
) const
{
    return
      - (
            this->alpha().boundaryField()[patchi]
           *this->thermo().kappa().boundaryField()[patchi]
           *this->thermo().T().boundaryField()[patchi].snGrad()
        );
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix>
Fourier<BasicThermophysicalTransportModel>::divq(volScalarField& he) const
{
    const thermoModel& thermo = this->thermo();

    // Return heat flux source as an implicit energy correction
    // to the temperature gradient flux
    return
        -fvc::laplacian(this->alpha()*thermo.kappa(), thermo.T())
        -fvm::laplacianCorrection
         (
             this->alpha()*thermo.kappa()/thermo.Cpv(),
             he
         );
}


template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField> Fourier<BasicThermophysicalTransportModel>::j
(
    const volScalarField& Yi
) const
{
    FatalErrorInFunction
        << type() << " supports single component systems only, " << nl
        << "    for multi-component transport select"
           " unityLewisFourier"
        << exit(FatalError);

    return tmp<surfaceScalarField>(nullptr);
}


template<class BasicThermophysicalTransportModel>
tmp<scalarField> Fourier<BasicThermophysicalTransportModel>::j
(
    const volScalarField& Yi,
    const label patchi
) const
{
    FatalErrorInFunction
        << type() << " supports single component systems only, " << nl
        << "    for multi-component transport select"
           " unityLewisFourier"
        << exit(FatalError);

    return tmp<scalarField>(nullptr);
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix>
Fourier<BasicThermophysicalTransportModel>::divj(volScalarField& Yi) const
{
    FatalErrorInFunction
        << type() << " supports single component systems only, " << nl
        << "    for multi-component transport select"
           " unityLewisFourier"
        << exit(FatalError);

    return tmp<fvScalarMatrix>(nullptr);
}


template<class BasicThermophysicalTransportModel>
void Fourier<BasicThermophysicalTransportModel>::predict()
{
    laminarThermophysicalTransportModel
    <
        BasicThermophysicalTransportModel
    >::predict();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
