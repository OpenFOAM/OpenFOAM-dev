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

#include "compressibleInterPhaseThermophysicalTransportModel.H"
#include "compressibleInterPhaseTransportModel.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleInterPhaseThermophysicalTransportModel::
compressibleInterPhaseThermophysicalTransportModel
(
    const compressibleInterPhaseTransportModel& momentumTransport
)
:
    thermophysicalTransportModel(momentumTransport.mixture_.mesh(), word::null),
    momentumTransport_(momentumTransport)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::compressibleInterPhaseThermophysicalTransportModel::read()
{
    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::compressibleInterPhaseThermophysicalTransportModel::kappaEff() const
{
    const compressibleTwoPhaseVoFMixture& mixture_ =
        momentumTransport_.mixture_;

    if (momentumTransport_.twoPhaseTransport_)
    {
        return
            mixture_.alpha1()
           *(
                mixture_.thermo1().kappa()
              + mixture_.thermo1().rho()*mixture_.thermo1().Cp()
               *momentumTransport_.momentumTransport1_->nut()
            )
          + mixture_.alpha2()
           *(
                mixture_.thermo2().kappa()
              + mixture_.thermo2().rho()*mixture_.thermo2().Cp()
               *momentumTransport_.momentumTransport2_->nut()
            );
    }
    else
    {
        return
            mixture_.alpha1()
           *(
                mixture_.thermo1().kappa()
              + mixture_.thermo1().rho()*mixture_.thermo1().Cp()
               *momentumTransport_.mixtureMomentumTransport_->nut()
            )
          + mixture_.alpha2()
           *(
                mixture_.thermo2().kappa()
              + mixture_.thermo2().rho()*mixture_.thermo2().Cp()
               *momentumTransport_.mixtureMomentumTransport_->nut()
            );
    }
}


Foam::tmp<Foam::scalarField>
Foam::compressibleInterPhaseThermophysicalTransportModel::kappaEff
(
    const label patchi
) const
{
    const compressibleTwoPhaseVoFMixture& mixture_ =
        momentumTransport_.mixture_;

    if (momentumTransport_.twoPhaseTransport_)
    {
        return
            mixture_.alpha1().boundaryField()[patchi]
           *(
                mixture_.thermo1().kappa().boundaryField()[patchi]
              + mixture_.thermo1().rho(patchi)
               *mixture_.thermo1().Cp().boundaryField()[patchi]
               *momentumTransport_.momentumTransport1_->nut(patchi)
            )
          + mixture_.alpha2().boundaryField()[patchi]
           *(
                mixture_.thermo2().kappa().boundaryField()[patchi]
              + mixture_.thermo2().rho(patchi)
               *mixture_.thermo2().Cp().boundaryField()[patchi]
               *momentumTransport_.momentumTransport2_->nut(patchi)
            );
    }
    else
    {
        return
            mixture_.alpha1().boundaryField()[patchi]
           *(
                mixture_.thermo1().kappa().boundaryField()[patchi]
              + mixture_.thermo1().rho(patchi)
               *mixture_.thermo1().Cp().boundaryField()[patchi]
               *momentumTransport_.mixtureMomentumTransport_->nut(patchi)
            )
          + mixture_.alpha2().boundaryField()[patchi]
           *(
                mixture_.thermo2().kappa().boundaryField()[patchi]
              + mixture_.thermo2().rho(patchi)
               *mixture_.thermo2().Cp().boundaryField()[patchi]
               *momentumTransport_.mixtureMomentumTransport_->nut(patchi)
            );
    }
}


Foam::tmp<Foam::volScalarField>
Foam::compressibleInterPhaseThermophysicalTransportModel::alphaEff() const
{
    const compressibleTwoPhaseVoFMixture& mixture_ =
        momentumTransport_.mixture_;

    if (momentumTransport_.twoPhaseTransport_)
    {
        return
            mixture_.alpha1()
           *(
                mixture_.thermo1().kappa()
              + mixture_.thermo1().rho()*mixture_.thermo1().Cp()
               *momentumTransport_.momentumTransport1_->nut()
            )/mixture_.thermo1().Cv()
          + mixture_.alpha2()
           *(
                mixture_.thermo2().kappa()
              + mixture_.thermo2().rho()*mixture_.thermo2().Cp()
               *momentumTransport_.momentumTransport2_->nut()
            )/mixture_.thermo2().Cv();
    }
    else
    {
        return
            mixture_.alpha1()
           *(
                mixture_.thermo1().kappa()
              + mixture_.thermo1().rho()*mixture_.thermo1().Cp()
               *momentumTransport_.mixtureMomentumTransport_->nut()
            )/mixture_.thermo1().Cv()
          + mixture_.alpha2()
           *(
                mixture_.thermo2().kappa()
              + mixture_.thermo2().rho()*mixture_.thermo2().Cp()
               *momentumTransport_.mixtureMomentumTransport_->nut()
            )/mixture_.thermo2().Cv();
    }
}


Foam::tmp<Foam::surfaceScalarField>
Foam::compressibleInterPhaseThermophysicalTransportModel::q() const
{
    return surfaceScalarField::New
    (
        "q",
        -fvc::interpolate(kappaEff())
        *fvc::snGrad(momentumTransport_.mixture_.T())
    );
}


Foam::tmp<Foam::scalarField>
Foam::compressibleInterPhaseThermophysicalTransportModel::qCorr
(
    const label patchi
) const
{
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::compressibleInterPhaseThermophysicalTransportModel::divq
(
    volScalarField& he
) const
{
    NotImplemented;

    return tmp<fvScalarMatrix>(nullptr);
}


void Foam::compressibleInterPhaseThermophysicalTransportModel::predict()
{}


void Foam::compressibleInterPhaseThermophysicalTransportModel::correct()
{}


// ************************************************************************* //
