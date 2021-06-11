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

#include "nonUnityLewisEddyDiffusivity.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvmSup.H"
#include "surfaceInterpolate.H"

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
    unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
    (
        typeName,
        momentumTransport,
        thermo,
        false
    ),

    Sct_("Sct", dimless, this->coeffDict_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
bool
nonUnityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::read()
{
    if
    (
        unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
        ::read()
    )
    {
        Sct_.read(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class TurbulenceThermophysicalTransportModel>
tmp<surfaceScalarField>
nonUnityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::q() const
{
    tmp<surfaceScalarField> tmpq
    (
        surfaceScalarField::New
        (
            IOobject::groupName
            (
                "q",
                this->momentumTransport().alphaRhoPhi().group()
            ),
           -fvc::interpolate(this->alpha()*this->kappaEff())
           *fvc::snGrad(this->thermo().T())
        )
    );

    const basicSpecieMixture& composition = this->thermo().composition();
    const PtrList<volScalarField>& Y = composition.Y();

    if (Y.size())
    {
        surfaceScalarField hGradY
        (
            surfaceScalarField::New
            (
                "hGradY",
                Y[0].mesh(),
                dimensionedScalar(dimEnergy/dimMass/dimLength, 0)
            )
        );

        forAll(Y, i)
        {
            const volScalarField hi
            (
                composition.Hs(i, this->thermo().p(), this->thermo().T())
            );

            hGradY += fvc::interpolate(hi)*fvc::snGrad(Y[i]);
        }

        tmpq.ref() -=
            fvc::interpolate
            (
                this->alpha()
               *this->thermo().alphaEff((this->Prt_/Sct_)*this->alphat())
            )*hGradY;
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
    tmp<fvScalarMatrix> tmpDivq
    (
        fvm::Su
        (
            -fvc::laplacian(this->alpha()*this->kappaEff(), this->thermo().T()),
            he
        )
    );

    const basicSpecieMixture& composition = this->thermo().composition();
    const PtrList<volScalarField>& Y = composition.Y();

    tmpDivq.ref() -=
        correction(fvm::laplacian(this->alpha()*this->alphaEff(), he));

    surfaceScalarField hGradY
    (
        surfaceScalarField::New
        (
            "hGradY",
            he.mesh(),
            dimensionedScalar(he.dimensions()/dimLength, 0)
        )
    );

    forAll(Y, i)
    {
        const volScalarField hi
        (
            composition.Hs(i, this->thermo().p(), this->thermo().T())
        );

        hGradY += fvc::interpolate(hi)*fvc::snGrad(Y[i]);
    }

    tmpDivq.ref() -=
        fvc::div
        (
            fvc::interpolate
            (
                this->alpha()
               *this->thermo().alphaEff((this->Prt_/Sct_)*this->alphat())
            )*hGradY*he.mesh().magSf()
        );

    return tmpDivq;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
