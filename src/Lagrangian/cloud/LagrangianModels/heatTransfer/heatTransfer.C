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

#include "heatTransfer.H"
#include "addToRunTimeSelectionTable.H"
#include "coupledToThermalFluid.H"
#include "LagrangianmSp.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
    defineTypeNameAndDebug(heatTransfer, 0);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Lagrangian::heatTransfer::heatTransfer
(
    const word& name,
    const LagrangianMesh& mesh,
    const dictionary& modelDict,
    const dictionary& stateDict
)
:
    LagrangianModel(name, mesh),
    cloudLagrangianModel(static_cast<const LagrangianModel&>(*this)),
    thermalCloud_(cloud<clouds::thermal>()),
    H
    (
        cloud().derivedField<scalar>
        (
            "H",
            *this,
            &heatTransfer::calcH
        )
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::wordList Foam::Lagrangian::heatTransfer::addSupFields() const
{
    return wordList({thermalCloud_.e.name()});
}


bool Foam::Lagrangian::heatTransfer::addsSupToField
(
    const word& fieldName,
    const word& eqnFieldName
) const
{
    return
        fieldName == thermalCloud_.e.name()
     && (
            eqnFieldName == thermalCloud_.e.name()
         || eqnFieldName == cloud<clouds::coupledToThermalFluid>().hec.name()
        );
}


void Foam::Lagrangian::heatTransfer::addSup
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& vOrM,
    const LagrangianSubScalarSubField& e,
    LagrangianEqn<scalar>& eqn
) const
{
    const clouds::coupledToThermalFluid& cctfCloud =
        cloud<clouds::coupledToThermalFluid>();

    const LagrangianSubMesh& subMesh = deltaT.mesh();

    const LagrangianSubScalarField& H = this->H(subMesh);

    const LagrangianSubScalarSubField T(subMesh.sub(thermalCloud_.T));
    const LagrangianSubScalarField& Tc = cctfCloud.Tc(subMesh);

    if (eqn.isPsi(e))
    {
        const LagrangianSubScalarSubField Cv(subMesh.sub(thermalCloud_.Cv));
        const LagrangianSubScalarField HbyCv(H/Cv);

        eqn.Su += H*(Tc - T) + HbyCv*e;
        eqn -= Lagrangianm::Sp(HbyCv, e);
    }
    else
    {
        const LagrangianSubScalarField& hec = cctfCloud.hec(subMesh);
        const LagrangianSubScalarField& Cpvc = cctfCloud.Cpvc(subMesh);
        const LagrangianSubScalarField HbyCpvc(H/Cpvc);

        eqn.Su += H*(Tc - T) - HbyCpvc*hec;
        eqn += Lagrangianm::Sp(HbyCpvc, hec);
    }
}


// ************************************************************************* //
