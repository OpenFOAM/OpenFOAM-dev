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

#include "RASeddyDiffusivity.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASThermophysicalTransportModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
void eddyDiffusivity<BasicThermophysicalTransportModel>::correctAlphat()
{
    alphat_ =
        this->momentumTransport().rho()
       *this->momentumTransport().nut()/Prt_;
    alphat_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
eddyDiffusivity<BasicThermophysicalTransportModel>::eddyDiffusivity
(
    const momentumTransportModel& momentumTransport
)
:
    RASThermophysicalTransportModel<BasicThermophysicalTransportModel>
    (
        typeName,
        momentumTransport
    ),

    Prt_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Prt",
            this->coeffDict_,
            1
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

template<class BasicThermophysicalTransportModel>
bool eddyDiffusivity<BasicThermophysicalTransportModel>::read()
{
    if
    (
        RASThermophysicalTransportModel
        <
            BasicThermophysicalTransportModel
        >::read()
    )
    {
        Prt_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicThermophysicalTransportModel>
tmp<volVectorField>eddyDiffusivity<BasicThermophysicalTransportModel>::q() const
{
    return volVectorField::New
    (
        IOobject::groupName
        (
            "q",
            this->momentumTransport().alphaRhoPhi().group()
        ),
       -this->alphaEff()
       *fvc::grad(this->thermo().he())
    );
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix>
eddyDiffusivity<BasicThermophysicalTransportModel>::divq
(
    volScalarField& he
) const
{
    return -fvm::laplacian(this->alphaEff(), he);
}


template<class BasicThermophysicalTransportModel>
void eddyDiffusivity<BasicThermophysicalTransportModel>::correct()
{
    RASThermophysicalTransportModel
    <
        BasicThermophysicalTransportModel
    >::correct();
    correctAlphat();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
