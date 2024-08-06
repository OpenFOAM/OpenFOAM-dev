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

#include "unityLewisEddyDiffusivity.H"
#include "fvcSnGrad.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceThermophysicalTransportModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
void unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::
correctAlphat()
{
    alphat_ =
        this->momentumTransport().rho()
       *this->momentumTransport().nut()/Prt_;
    alphat_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::
unityLewisEddyDiffusivity
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo,
    const bool allowDefaultPrt
)
:
    unityLewisEddyDiffusivity
    (
        typeName,
        momentumTransport,
        thermo,
        allowDefaultPrt
    )
{}


template<class TurbulenceThermophysicalTransportModel>
unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::
unityLewisEddyDiffusivity
(
    const word& type,
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo,
    const bool allowDefaultPrt
)
:
    TurbulenceThermophysicalTransportModel
    (
        type,
        momentumTransport,
        thermo
    ),

    Prt_
    (
        allowDefaultPrt
      ? dimensioned<scalar>("Prt", dimless, this->coeffDict(), 1)
      : dimensioned<scalar>("Prt", dimless, this->coeffDict())
    ),

    alphat_
    (
        IOobject
        (
            IOobject::groupName
            (
                "alphat",
                this->momentumTransport().alphaRhoPhi().group()
            ),
            momentumTransport.time().name(),
            momentumTransport.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        momentumTransport.mesh()
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
bool unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::read()
{
    if (TurbulenceThermophysicalTransportModel::read())
    {
        Prt_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class TurbulenceThermophysicalTransportModel>
tmp<surfaceScalarField>
unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::q() const
{
    return surfaceScalarField::New
    (
        IOobject::groupName
        (
            "q",
            this->momentumTransport().alphaRhoPhi().group()
        ),
       -fvc::interpolate(this->alphaEff()*this->alpha())
       *fvc::snGrad(this->thermo().he())
    );
}


template<class TurbulenceThermophysicalTransportModel>
tmp<fvScalarMatrix>
unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::divq
(
    volScalarField& he
) const
{
    return -fvm::laplacian(this->alpha()*this->alphaEff(), he);
}


template<class TurbulenceThermophysicalTransportModel>
tmp<surfaceScalarField>
unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::j
(
    const volScalarField& Yi
) const
{
    return surfaceScalarField::New
    (
        IOobject::groupName
        (
            "j(" + Yi.name() + ')',
            this->momentumTransport().alphaRhoPhi().group()
        ),
       -fvc::interpolate(this->DEff(Yi)*this->alpha())*fvc::snGrad(Yi)
    );
}


template<class TurbulenceThermophysicalTransportModel>
tmp<fvScalarMatrix>
unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::divj
(
    volScalarField& Yi
) const
{
    return -fvm::laplacian(this->alpha()*this->DEff(Yi), Yi);
}


template<class TurbulenceThermophysicalTransportModel>
void unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::
predict()
{
    TurbulenceThermophysicalTransportModel::predict();
    correctAlphat();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
