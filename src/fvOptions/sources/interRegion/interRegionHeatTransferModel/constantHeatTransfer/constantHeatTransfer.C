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

#include "constantHeatTransfer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(constantHeatTransfer, 0);
    addToRunTimeSelectionTable
    (
        option,
        constantHeatTransfer,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::constantHeatTransfer::constantHeatTransfer
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    interRegionHeatTransferModel(name, modelType, dict, mesh),
    htcConst_(),
    AoV_()
{
    if (active() && master_)
    {
        htcConst_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "htcConst",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );

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

        htc_ = htcConst_()*AoV_();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::constantHeatTransfer::~constantHeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::constantHeatTransfer::calculateHtc()
{
    // do nothing
}


void Foam::fv::constantHeatTransfer::writeData(Ostream& os) const
{
    os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
    interRegionHeatTransferModel::writeData(os);

    os << indent << type() + "Coeffs" << nl;

    coeffs_.write(os);

    os << decrIndent << indent << token::END_BLOCK << endl;
}


bool Foam::fv::constantHeatTransfer::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
