/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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

#include "MulticomponentPhaseModel.H"

#include "phaseSystem.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MulticomponentPhaseModel<BasePhaseModel>::MulticomponentPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index)
{
    PtrList<volScalarField>& Y = this->thermo_->composition().Y();

    forAll(Y, i)
    {
        if (this->thermo_->composition().solve(i))
        {
            const label j = YActive_.size();
            YActive_.resize(j + 1);
            YActive_.set(j, &Y[i]);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MulticomponentPhaseModel<BasePhaseModel>::~MulticomponentPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::MulticomponentPhaseModel<BasePhaseModel>::correctSpecies()
{
    this->thermo_->composition().normalise();
    BasePhaseModel::correctSpecies();
}


template<class BasePhaseModel>
bool Foam::MulticomponentPhaseModel<BasePhaseModel>::pure() const
{
    return false;
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MulticomponentPhaseModel<BasePhaseModel>::YiEqn(volScalarField& Yi)
{
    const volScalarField& alpha = *this;
    const volScalarField& rho = this->rho();

    const tmp<surfaceScalarField> talphaRhoPhi(this->alphaRhoPhi());
    const surfaceScalarField& alphaRhoPhi(talphaRhoPhi());

    return
    (
        fvm::ddt(alpha, rho, Yi)
      + fvm::div(alphaRhoPhi, Yi, "div(" + alphaRhoPhi.name() + ",Yi)")
      + this->divj(Yi)
     ==
        alpha*this->R(Yi)

      - correction
        (
            fvm::Sp
            (
                max(this->residualAlpha() - alpha, scalar(0))*rho
               /this->mesh().time().deltaT(),
                Yi
            )
        )
    );
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::MulticomponentPhaseModel<BasePhaseModel>::Y() const
{
    return this->thermo_->composition().Y();
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::MulticomponentPhaseModel<BasePhaseModel>::Y(const word& name) const
{
    return this->thermo_->composition().Y(name);
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::MulticomponentPhaseModel<BasePhaseModel>::YRef()
{
    return this->thermo_->composition().Y();
}


template<class BasePhaseModel>
const Foam::UPtrList<Foam::volScalarField>&
Foam::MulticomponentPhaseModel<BasePhaseModel>::YActive() const
{
    return YActive_;
}


template<class BasePhaseModel>
Foam::UPtrList<Foam::volScalarField>&
Foam::MulticomponentPhaseModel<BasePhaseModel>::YActiveRef()
{
    return YActive_;
}


// ************************************************************************* //
