/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    eddyDiffusivity
    (
        typeName,
        momentumTransport,
        thermo,
        false
    )
{}


template<class TurbulenceThermophysicalTransportModel>
eddyDiffusivity<TurbulenceThermophysicalTransportModel>::eddyDiffusivity
(
    const word& type,
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo,
    const bool allowDefaultPrt
)
:
    TurbulenceThermophysicalTransportModel
    (
        typeName,
        momentumTransport,
        thermo
    ),

    Prt_
    (
        allowDefaultPrt
      ? dimensioned<scalar>::lookupOrAddToDict
        (
            "Prt",
            this->coeffDict_,
            1
        )
      : dimensioned<scalar>
        (
            "Prt",
            dimless,
            this->coeffDict_
        )
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
        Prt_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class TurbulenceThermophysicalTransportModel>
tmp<volVectorField>
eddyDiffusivity<TurbulenceThermophysicalTransportModel>::q() const
{
    return volVectorField::New
    (
        IOobject::groupName
        (
            "q",
            this->momentumTransport().alphaRhoPhi().group()
        ),
       -this->alphaEff()*this->alpha()*fvc::grad(this->thermo().he())
    );
}


template<class TurbulenceThermophysicalTransportModel>
tmp<fvScalarMatrix>
eddyDiffusivity<TurbulenceThermophysicalTransportModel>::divq
(
    volScalarField& he
) const
{
    return -fvm::laplacian(this->alpha()*this->alphaEff(), he);
}


template<class TurbulenceThermophysicalTransportModel>
tmp<volVectorField>
eddyDiffusivity<TurbulenceThermophysicalTransportModel>::j
(
    const volScalarField& Yi
) const
{
    return volVectorField::New
    (
        IOobject::groupName
        (
            "j(" + Yi.name() + ')',
            this->momentumTransport().alphaRhoPhi().group()
        ),
       -this->DEff(Yi)*this->alpha()*fvc::grad(Yi)
    );
}


template<class TurbulenceThermophysicalTransportModel>
tmp<fvScalarMatrix>
eddyDiffusivity<TurbulenceThermophysicalTransportModel>::divj
(
    volScalarField& Yi
) const
{
    return -fvm::laplacian(this->alpha()*this->DEff(Yi), Yi);
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
