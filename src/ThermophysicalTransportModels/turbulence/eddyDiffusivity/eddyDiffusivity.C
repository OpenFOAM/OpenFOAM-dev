/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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

#include "eddyDiffusivity.H"
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceThermophysicalTransportModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
void eddyDiffusivity<TurbulenceThermophysicalTransportModel>::correctAlphat()
{
    alphat_ =
        this->momentumTransport().rho()
       *this->momentumTransport().nut()/Prt_;
    alphat_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
eddyDiffusivity<TurbulenceThermophysicalTransportModel>::eddyDiffusivity
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    TurbulenceThermophysicalTransportModel
    (
        typeName,
        momentumTransport,
        thermo
    ),

    Prt_("Prt", dimless, this->coeffDict_),

    alphat_
    (
        IOobject
        (
            IOobject::groupName
            (
                "alphat",
                this->momentumTransport().alphaRhoPhi().group()
            ),
            momentumTransport.time().timeName(),
            momentumTransport.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        momentumTransport.mesh()
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
bool eddyDiffusivity<TurbulenceThermophysicalTransportModel>::read()
{
    if (TurbulenceThermophysicalTransportModel::read())
    {
        Prt_.read(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class TurbulenceThermophysicalTransportModel>
tmp<volScalarField>
eddyDiffusivity<TurbulenceThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi
) const
{
    FatalErrorInFunction
        << type() << " supports single component systems only, " << nl
        << "    for multi-component transport select"
           " nonUnityLewisEddyDiffusivity or unityLewisEddyDiffusivity"
        << exit(FatalError);

    return tmp<volScalarField>(nullptr);
}


template<class TurbulenceThermophysicalTransportModel>
tmp<scalarField>
eddyDiffusivity<TurbulenceThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{
    FatalErrorInFunction
        << type() << " supports single component systems only, " << nl
        << "    for multi-component transport select"
           " nonUnityLewisEddyDiffusivity or unityLewisEddyDiffusivity"
        << exit(FatalError);

    return tmp<scalarField>(nullptr);
}


template<class TurbulenceThermophysicalTransportModel>
tmp<surfaceScalarField>
eddyDiffusivity<TurbulenceThermophysicalTransportModel>::q() const
{
    return surfaceScalarField::New
    (
        IOobject::groupName
        (
            "q",
            this->momentumTransport().alphaRhoPhi().group()
        ),
       -fvc::interpolate(this->alpha()*this->kappaEff())
       *fvc::snGrad(this->thermo().T())
    );
}


template<class TurbulenceThermophysicalTransportModel>
tmp<fvScalarMatrix>
eddyDiffusivity<TurbulenceThermophysicalTransportModel>::divq
(
    volScalarField& he
) const
{
    // Return heat flux source as an implicit energy correction
    // to the temperature gradient flux
    return
        -correction(fvm::laplacian(this->alpha()*this->alphaEff(), he))
        -fvc::laplacian(this->alpha()*this->kappaEff(), this->thermo().T());
}


template<class TurbulenceThermophysicalTransportModel>
tmp<surfaceScalarField>
eddyDiffusivity<TurbulenceThermophysicalTransportModel>::j
(
    const volScalarField& Yi
) const
{
    FatalErrorInFunction
        << type() << " supports single component systems only, " << nl
        << "    for multi-component transport select"
           " nonUnityLewisEddyDiffusivity or unityLewisEddyDiffusivity"
        << exit(FatalError);

    return tmp<surfaceScalarField>(nullptr);
}


template<class TurbulenceThermophysicalTransportModel>
tmp<fvScalarMatrix>
eddyDiffusivity<TurbulenceThermophysicalTransportModel>::divj
(
    volScalarField& Yi
) const
{
    FatalErrorInFunction
        << type() << " supports single component systems only, " << nl
        << "    for multi-component transport select"
           " nonUnityLewisEddyDiffusivity or unityLewisEddyDiffusivity"
        << exit(FatalError);

    return tmp<fvScalarMatrix>(nullptr);
}


template<class TurbulenceThermophysicalTransportModel>
void eddyDiffusivity<TurbulenceThermophysicalTransportModel>::correct()
{
    TurbulenceThermophysicalTransportModel::correct();
    correctAlphat();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
