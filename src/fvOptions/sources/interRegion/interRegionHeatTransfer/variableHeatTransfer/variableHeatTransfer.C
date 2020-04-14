/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "variableHeatTransfer.H"
#include "thermophysicalTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(variableHeatTransfer, 0);
    addToRunTimeSelectionTable
    (
        option,
        variableHeatTransfer,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::variableHeatTransfer::variableHeatTransfer
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    interRegionHeatTransferModel(name, modelType, dict, mesh),
    UNbrName_(coeffs_.lookupOrDefault<word>("UNbr", "U")),
    a_(0),
    b_(0),
    c_(0),
    ds_(0),
    Pr_(0),
    AoV_()
{
    if (master_)
    {
        a_ = coeffs_.lookup<scalar>("a");
        b_ = coeffs_.lookup<scalar>("b");
        c_ = coeffs_.lookup<scalar>("c");
        ds_ = coeffs_.lookup<scalar>("ds");
        Pr_ = coeffs_.lookup<scalar>("Pr");
        AoV_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "AoV",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::variableHeatTransfer::~variableHeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::variableHeatTransfer::calculateHtc()
{
    const fvMesh& nbrMesh =
        mesh_.time().lookupObject<fvMesh>(nbrRegionName());

    const thermophysicalTransportModel& nbrTtm =
        nbrMesh.lookupObject<thermophysicalTransportModel>
        (
            thermophysicalTransportModel::typeName
        );

    const compressibleMomentumTransportModel& nbrTurb =
        nbrTtm.momentumTransport();

    const fluidThermo& nbrThermo = nbrTtm.thermo();

    const volVectorField& UNbr =
        nbrMesh.lookupObject<volVectorField>(UNbrName_);

    const volScalarField ReNbr(mag(UNbr)*ds_*nbrThermo.rho()/nbrTurb.mut());

    const volScalarField NuNbr(a_*pow(ReNbr, b_)*pow(Pr_, c_));

    const scalarField htcNbr(NuNbr*nbrTtm.kappaEff()/ds_);

    const scalarField htcNbrMapped(interpolate(htcNbr));

    htc_.primitiveFieldRef() = htcNbrMapped*AoV_;
}


bool Foam::fv::variableHeatTransfer::read(const dictionary& dict)
{
    if (interRegionHeatTransferModel::read(dict))
    {
        coeffs_.readIfPresent("UNbr", UNbrName_);

        coeffs_.readIfPresent("a", a_);
        coeffs_.readIfPresent("b", b_);
        coeffs_.readIfPresent("c", c_);
        coeffs_.readIfPresent("ds", ds_);
        coeffs_.readIfPresent("Pr", Pr_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
