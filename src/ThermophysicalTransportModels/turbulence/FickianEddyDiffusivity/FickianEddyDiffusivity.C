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

#include "FickianEddyDiffusivity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceThermophysicalTransportModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
FickianEddyDiffusivity<TurbulenceThermophysicalTransportModel>::
FickianEddyDiffusivity
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    Fickian<unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>>
    (
        typeName,
        momentumTransport,
        thermo
    ),

    Sct_("Sct", dimless, this->coeffDict_)
{
    read();
    this->correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
bool
FickianEddyDiffusivity<TurbulenceThermophysicalTransportModel>::read()
{
    if
    (
        Fickian
        <
            unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
        >::read()
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
tmp<volScalarField>
FickianEddyDiffusivity<TurbulenceThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi
) const
{
    return volScalarField::New
    (
        "DEff",
        Fickian
        <
            unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
        >::DEff(Yi)
      + (this->Prt_/Sct_)*this->alphat()
    );
}


template<class TurbulenceThermophysicalTransportModel>
tmp<scalarField>
FickianEddyDiffusivity<TurbulenceThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{
    return
        Fickian
        <
            unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
        >::DEff(Yi, patchi)
      + this->Prt_.value()/Sct_.value()*this->alphat(patchi);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
