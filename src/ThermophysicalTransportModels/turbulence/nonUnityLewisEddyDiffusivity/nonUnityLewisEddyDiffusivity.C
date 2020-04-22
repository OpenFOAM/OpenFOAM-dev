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

#include "nonUnityLewisEddyDiffusivity.H"
#include "fvcLaplacian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceThermophysicalTransportModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
nonUnityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::
nonUnityLewisEddyDiffusivity
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    eddyDiffusivity<TurbulenceThermophysicalTransportModel>
    (
        typeName,
        momentumTransport,
        thermo,
        false
    ),

    Sct_
    (
        dimensioned<scalar>
        (
            "Sct",
            dimless,
            this->coeffDict_
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
bool
nonUnityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::read()
{
    if (eddyDiffusivity<TurbulenceThermophysicalTransportModel>::read())
    {
        Sct_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class TurbulenceThermophysicalTransportModel>
tmp<volVectorField>
nonUnityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::q() const
{
    tmp<volVectorField> tmpq =
        eddyDiffusivity<TurbulenceThermophysicalTransportModel>::q();

    if (mag(this->Prt_ - Sct_).value() > small)
    {
        const basicSpecieMixture& composition = this->thermo().composition();

        const PtrList<volScalarField>& Y = composition.Y();

        volScalarField alphatMinusDt
        (
            "alphatMinusDt",
            this->alphat()*(1 - this->Prt_/Sct_)
        );

        forAll(Y, i)
        {
            tmpq.ref() +=
                this->alpha()
               *alphatMinusDt
               *composition.HE(i, this->thermo().p(), this->thermo().T())
               *fvc::grad(Y[i]);
        }
    }

    return tmpq;
}


template<class TurbulenceThermophysicalTransportModel>
tmp<fvScalarMatrix>
nonUnityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::divq
(
    volScalarField& he
) const
{
    tmp<fvScalarMatrix> tmpDivq =
        eddyDiffusivity<TurbulenceThermophysicalTransportModel>::divq(he);

    if (mag(this->Prt_ - Sct_).value() > small)
    {
        const basicSpecieMixture& composition = this->thermo().composition();

        const PtrList<volScalarField>& Y = composition.Y();

        volScalarField alphatMinusDt
        (
            "alphatMinusDt",
            this->alphat()*(1 - this->Prt_/Sct_)
        );

        forAll(Y, i)
        {
            tmpDivq.ref() +=
                fvc::laplacian
                (
                    this->alpha()
                   *alphatMinusDt
                   *composition.HE(i, this->thermo().p(), this->thermo().T()),
                    Y[i]
                );
        }
    }

    return tmpDivq;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
