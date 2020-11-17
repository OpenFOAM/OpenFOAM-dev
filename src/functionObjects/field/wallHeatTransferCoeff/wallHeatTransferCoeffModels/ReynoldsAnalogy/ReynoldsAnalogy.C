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

#include "ReynoldsAnalogy.H"
#include "kinematicMomentumTransportModel.H"
#include "fluidThermoMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallHeatTransferCoeffModels
{
    defineTypeNameAndDebug(ReynoldsAnalogy, 0);
    addToRunTimeSelectionTable
    (
        wallHeatTransferCoeffModel,
        ReynoldsAnalogy,
        word
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallHeatTransferCoeffModels::ReynoldsAnalogy::ReynoldsAnalogy
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    wallHeatTransferCoeffModel(name, mesh, dict),
    mesh_(mesh),
    Uref_("Uref", dimVelocity, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallHeatTransferCoeffModels::ReynoldsAnalogy::~ReynoldsAnalogy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::wallHeatTransferCoeffModels::ReynoldsAnalogy::read
(
    const dictionary& dict
)
{
    Uref_.read(dict);
    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::wallHeatTransferCoeffModels::ReynoldsAnalogy::htcByRhoCp
(
    const momentumTransportModel& mmtm,
    const labelHashSet& patches
) const
{
    tmp<volSymmTensorField> ttau(this->tau(mmtm, mesh_));
    const volSymmTensorField::Boundary& tauBf = ttau.ref().boundaryField();

    // Create temporary field for heat transfer coefficient
    tmp<volScalarField> thtcByRhoCp
    (
        volScalarField::New
        (
            type(),
            mesh_,
            dimensionedScalar(dimLength/dimTime, 0)
        )
    );

    volScalarField::Boundary& thtcByRhoCpBf =
        thtcByRhoCp.ref().boundaryFieldRef();

    forAll(thtcByRhoCpBf, patchi)
    {
        if (!thtcByRhoCpBf[patchi].coupled())
        {
            const vectorField nf(tauBf[patchi].patch().nf());

            // Wall shear stress in [m^2/s^2]
            tmp<vectorField> tauwp(-nf&tauBf[patchi]);

            // Non-dimensional skin friction coefficient [-]
            const scalarField Cf(2*mag(tauwp)/sqr(Uref_.value()));

            thtcByRhoCpBf[patchi] = 0.5*Uref_.value()*Cf;
        }
    }

    return thtcByRhoCp;
}

// ************************************************************************* //
