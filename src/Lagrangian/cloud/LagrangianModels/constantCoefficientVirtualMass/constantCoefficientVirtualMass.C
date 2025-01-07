/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "constantCoefficientVirtualMass.H"
#include "addToRunTimeSelectionTable.H"
#include "coupledToIncompressibleFluid.H"
#include "coupledToFluid.H"
#include "shaped.H"
#include "LagrangianmDdt.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
    defineTypeNameAndDebug(constantCoefficientVirtualMass, 0);
    addToRunTimeSelectionTable
    (
        LagrangianModel,
        constantCoefficientVirtualMass,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarField>
Foam::Lagrangian::constantCoefficientVirtualMass::calcVOrMAdd
(
    const LagrangianModelRef& model,
    const LagrangianSubMesh& subMesh
) const
{
    const LagrangianSubScalarField& v = cloud<clouds::shaped>().v(subMesh);

    return
        isCloud<clouds::coupledToIncompressibleFluid>()
      ? Cvm_*v/cloud<clouds::coupledToIncompressibleFluid>().rhoByRhoc
      : Cvm_*v*cloud<clouds::coupledToFluid>().rhoc(subMesh);
}


void Foam::Lagrangian::constantCoefficientVirtualMass::addUSup
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubVectorSubField& U,
    LagrangianEqn<vector>& eqn
) const
{
    const LagrangianSubScalarSubField& vOrMAdd = this->vOrMAdd.ref(U.mesh());

    if (eqn.isPsi(U))
    {
        eqn -= Lagrangianm::ddt(deltaT, vOrMAdd, U);
    }
    else
    {
        eqn -= Lagrangianm::ddt0(deltaT, vOrMAdd, U);
    }

    eqn.Su += vOrMAdd*cloud<clouds::coupled>().DUDtc(U.mesh());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Lagrangian::constantCoefficientVirtualMass::constantCoefficientVirtualMass
(
    const word& name,
    const LagrangianMesh& mesh,
    const dictionary& modelDict,
    const dictionary& stateDict
)
:
    LagrangianModel(name, mesh),
    cloudLagrangianModel(static_cast<const LagrangianModel&>(*this)),
    Cvm_("Cvm", dimless, modelDict),
    vOrMAdd
    (
        cloud().derivedField<scalar>
        (
            *this,
            &constantCoefficientVirtualMass::calcVOrMAdd
        )
    )
{
    cloud<clouds::coupled>().DUDtc.psi();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::wordList
Foam::Lagrangian::constantCoefficientVirtualMass::addSupFields() const
{
    return wordList(1, cloud().U.name());
}


void Foam::Lagrangian::constantCoefficientVirtualMass::addSup
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubVectorSubField& U,
    LagrangianEqn<vector>& eqn
) const
{
    assertCloud<clouds::coupledToIncompressibleFluid>();

    addUSup(deltaT, U, eqn);
}


void Foam::Lagrangian::constantCoefficientVirtualMass::addSup
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& vOrM,
    const LagrangianSubVectorSubField& U,
    LagrangianEqn<vector>& eqn
) const
{
    assertCloud<clouds::coupledToIncompressibleFluid, clouds::coupledToFluid>();

    addUSup(deltaT, U, eqn);
}


// ************************************************************************* //
