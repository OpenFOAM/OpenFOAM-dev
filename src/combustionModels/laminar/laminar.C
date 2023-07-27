/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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

#include "laminar.H"
#include "fvmSup.H"
#include "localEulerDdtScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(laminar, 0);
    addToRunTimeSelectionTable(combustionModel, laminar, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::laminar::laminar
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
:
    combustionModel
    (
        modelType,
        thermo,
        turb,
        combustionProperties
    ),
    integrateReactionRate_
    (
        this->coeffs().lookupOrDefault("integrateReactionRate", true)
    ),
    maxIntegrationTime_
    (
        this->coeffs().lookupOrDefault("maxIntegrationTime", vGreat)
    ),
    outerCorrect_
    (
        this->coeffs().lookupOrDefault("outerCorrect", false)
    ),
    timeIndex_(-1),
    chemistryPtr_(basicChemistryModel::New(thermo))
{
    if (integrateReactionRate_)
    {
        Info<< "    using integrated reaction rate" << endl;
    }
    else
    {
        Info<< "    using instantaneous reaction rate" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::laminar::~laminar()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::laminar::correct()
{
    if
    (
        (integrateReactionRate_ || !outerCorrect_)
     && timeIndex_ == this->mesh().time().timeIndex()
    )
    {
        return;
    }

    if (integrateReactionRate_)
    {
        if (fv::localEulerDdt::enabled(this->mesh()))
        {
            const scalarField& rDeltaT =
                fv::localEulerDdt::localRDeltaT(this->mesh());

            chemistryPtr_->solve(min(1/rDeltaT, maxIntegrationTime_)());
        }
        else
        {
            const scalar deltaT = this->mesh().time().deltaTValue();

            chemistryPtr_->solve(min(deltaT, maxIntegrationTime_));
        }
    }
    else
    {
        chemistryPtr_->calculate();
    }

    timeIndex_ = this->mesh().time().timeIndex();
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::combustionModels::laminar::R(const label speciei) const
{
    return chemistryPtr_->RR()[speciei];
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::laminar::R(volScalarField& Y) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));
    fvScalarMatrix& Su = tSu.ref();

    const label specieI = this->thermo().species()[Y.member()];
    Su += chemistryPtr_->RR()[specieI];

    return tSu;
}


Foam::tmp<Foam::volScalarField> Foam::combustionModels::laminar::Qdot() const
{
    return chemistryPtr_->Qdot();
}


bool Foam::combustionModels::laminar::read()
{
    if (combustionModel::read())
    {
        integrateReactionRate_ =
            this->coeffs().lookupOrDefault("integrateReactionRate", true);
        maxIntegrationTime_ =
            this->coeffs().lookupOrDefault("maxIntegrationTime", vGreat);
        outerCorrect_ =
            this->coeffs().lookupOrDefault("outerCorrect", false);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
