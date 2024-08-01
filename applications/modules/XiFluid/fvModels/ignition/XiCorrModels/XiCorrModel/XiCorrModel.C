/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "XiCorrModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(XiCorrModel, 0);
    defineRunTimeSelectionTable(XiCorrModel, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::XiCorrModel::readCoeffs(const dictionary& dict)
{
    bMin_.readIfPresent(dict);
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiCorrModel::XiCorrModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    set_(mesh, dict),
    bMin_("bMin", dimless, 0.001)
{
    readCoeffs(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiCorrModel::~XiCorrModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::XiCorrModel::XiCorr
(
    volScalarField& Xi,
    const volScalarField& b,
    const volScalarField& mgb
) const
{
    const labelUList cells = set_.cells();
    const scalarField bCells(b, cells);

    scalar bMin = min(bCells);
    reduce(bMin, minOp<scalar>());

    if (bMin > bMin_.value())
    {
        const scalarField Vcells(b.mesh().V(), cells);

        // Calculate volume of ignition kernel
        const dimensionedScalar Vk
        (
            "Vk",
            dimVolume,
            gSum((1 - bCells)*Vcells)
        );

        if (Vk.value() > small)
        {
            // Calculate kernel area from its volume
            const dimensionedScalar Ak(this->Ak(Vk));

            const scalarField mgbCells(mgb, cells);

            // Calculate kernel area from b field
            const dimensionedScalar AkEst(gSum(mgbCells*Vcells));

            const scalar XiCorr = max(min((Ak/AkEst).value(), 10.0), 1.0);

            Info<< "XiCorr = " << XiCorr << ", bMin = " << bMin << endl;

            Xi *= XiCorr;
        }
    }
}


void Foam::XiCorrModel::topoChange
(
    const polyTopoChangeMap& map
)
{
    set_.topoChange(map);
}


void Foam::XiCorrModel::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void Foam::XiCorrModel::distribute
(
    const polyDistributionMap& map
)
{
    set_.distribute(map);
}


bool Foam::XiCorrModel::movePoints()
{
    set_.movePoints();
    return true;
}


bool Foam::XiCorrModel::read(const dictionary& dict)
{
    const dictionary& XiCorrDict
    (
        dict.subDict("XiCorr").optionalSubDict(type() + "Coeffs")
    );

    set_.read(XiCorrDict);
    return readCoeffs(XiCorrDict);
}


// ************************************************************************* //
