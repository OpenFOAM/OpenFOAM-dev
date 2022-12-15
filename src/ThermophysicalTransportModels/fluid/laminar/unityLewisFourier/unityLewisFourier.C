/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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

#include "unityLewisFourier.H"
#include "fvmLaplacian.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
unityLewisFourier<BasicThermophysicalTransportModel>::unityLewisFourier
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    unityLewisFourier
    (
        typeName,
        momentumTransport,
        thermo
    )
{}


template<class BasicThermophysicalTransportModel>
unityLewisFourier<BasicThermophysicalTransportModel>::unityLewisFourier
(
    const word& type,
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    laminarThermophysicalTransportModel<BasicThermophysicalTransportModel>
    (
        type,
        momentumTransport,
        thermo
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
const dictionary&
unityLewisFourier<BasicThermophysicalTransportModel>::coeffDict() const
{
    return dictionary::null;
}


template<class BasicThermophysicalTransportModel>
bool unityLewisFourier<BasicThermophysicalTransportModel>::read()
{
    return true;
}


template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField>
unityLewisFourier<BasicThermophysicalTransportModel>::q() const
{
    return surfaceScalarField::New
    (
        IOobject::groupName
        (
            "q",
            this->momentumTransport().alphaRhoPhi().group()
        ),
       -fvc::interpolate
        (
            this->alpha()*this->thermo().kappa()/this->thermo().Cpv()
        )
       *fvc::snGrad(this->thermo().he())
    );
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix>
unityLewisFourier<BasicThermophysicalTransportModel>::
divq(volScalarField& he) const
{
    volScalarField alphahe
    (
        volScalarField::New
        (
            "alphahe",
            this->thermo().kappa()/this->thermo().Cpv()
        )
    );

    return -fvm::laplacian(this->alpha()*alphahe, he);
}


template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField>unityLewisFourier<BasicThermophysicalTransportModel>::j
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
       -fvc::interpolate(this->alpha()*this->DEff(Yi))
       *fvc::snGrad(Yi)
    );
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix>
unityLewisFourier<BasicThermophysicalTransportModel>::
divj(volScalarField& Yi) const
{
    return -fvm::laplacian(this->alpha()*this->DEff(Yi), Yi);
}


template<class BasicThermophysicalTransportModel>
void unityLewisFourier<BasicThermophysicalTransportModel>::predict()
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
