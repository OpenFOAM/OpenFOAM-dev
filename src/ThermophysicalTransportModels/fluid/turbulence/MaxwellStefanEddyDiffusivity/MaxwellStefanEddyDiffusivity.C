/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "MaxwellStefanEddyDiffusivity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceThermophysicalTransportModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
MaxwellStefanEddyDiffusivity<TurbulenceThermophysicalTransportModel>::
MaxwellStefanEddyDiffusivity
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo,
    const word& type
)
:
    MaxwellStefan
    <
        unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
    >
    (
        type,
        momentumTransport,
        thermo
    ),

    Sct_("Sct", dimless, this->typeDict(type))
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
bool
MaxwellStefanEddyDiffusivity<TurbulenceThermophysicalTransportModel>::read()
{
    if
    (
        MaxwellStefan
        <
            unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
        >::read()
    )
    {
        Sct_.read(this->typeDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class TurbulenceThermophysicalTransportModel>
tmp<volScalarField>
MaxwellStefanEddyDiffusivity<TurbulenceThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi
) const
{
    return volScalarField::New
    (
        "DEff",
        MaxwellStefan
        <
            unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
        >::DEff(Yi)
      + (this->Prt_/Sct_)*this->alphat()
    );
}


template<class TurbulenceThermophysicalTransportModel>
tmp<scalarField>
MaxwellStefanEddyDiffusivity<TurbulenceThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{
    return
        MaxwellStefan
        <
            unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
        >::DEff(Yi, patchi)
      + this->Prt_.value()/Sct_.value()*this->alphat(patchi);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
