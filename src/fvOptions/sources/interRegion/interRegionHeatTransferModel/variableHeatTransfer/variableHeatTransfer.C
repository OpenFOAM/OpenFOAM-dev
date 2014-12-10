/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
#include "turbulenceModel.H"
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
    UNbrName_(coeffs_.lookupOrDefault<word>("UNbrName", "U")),
    a_(0),
    b_(0),
    c_(0),
    ds_(0),
    Pr_(0),
    AoV_()
{
    if (master_)
    {
        a_ = readScalar(coeffs_.lookup("a"));
        b_ = readScalar(coeffs_.lookup("b"));
        c_ = readScalar(coeffs_.lookup("c"));
        ds_ = readScalar(coeffs_.lookup("ds"));
        Pr_ = readScalar(coeffs_.lookup("Pr"));
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

    const compressible::turbulenceModel& nbrTurb =
        nbrMesh.lookupObject<compressible::turbulenceModel>("turbulenceModel");

    const fluidThermo& nbrThermo =
        nbrMesh.lookupObject<fluidThermo>("thermophysicalProperties");

    const volVectorField& UNbr =
        nbrMesh.lookupObject<volVectorField>(UNbrName_);

    const volScalarField ReNbr(mag(UNbr)*ds_*nbrThermo.rho()/nbrTurb.mut());

    const volScalarField NuNbr(a_*pow(ReNbr, b_)*pow(Pr_, c_));

    const scalarField htcNbr(NuNbr*nbrTurb.kappaEff()/ds_);

    const scalarField htcNbrMapped(interpolate(htcNbr));

    htc_.internalField() = htcNbrMapped*AoV_;
}


void Foam::fv::variableHeatTransfer::writeData(Ostream& os) const
{
    os  << indent << token::BEGIN_BLOCK << incrIndent << nl;

    interRegionHeatTransferModel::writeData(os);

    os << indent << type() + "Coeffs" << nl;

    coeffs_.write(os);

    os << decrIndent << indent << token::END_BLOCK << endl;
}


bool Foam::fv::variableHeatTransfer::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        coeffs_.readIfPresent("UNbrName", UNbrName_);

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
