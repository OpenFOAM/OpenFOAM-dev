/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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

#include "waveAtmBoundaryLayerSuperposition.H"
#include "uniformDimensionedFields.H"
#include "atmBoundaryLayer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(waveAtmBoundaryLayerSuperposition, 0);
    addToRunTimeSelectionTable
    (
        waveSuperposition,
        waveAtmBoundaryLayerSuperposition,
        objectRegistry
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveAtmBoundaryLayerSuperposition::waveAtmBoundaryLayerSuperposition
(
    const objectRegistry& db
)
:
    waveSuperposition(db),
    UGasRef_(lookup("UGasRef")),
    hRef_(lookup<scalar>("hRef")),
    hWaveMin_(lookup<scalar>("hWaveMin")),
    hWaveMax_(lookup<scalar>("hWaveMax"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveAtmBoundaryLayerSuperposition::~waveAtmBoundaryLayerSuperposition()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::waveAtmBoundaryLayerSuperposition::UGas
(
    const scalar t,
    const vectorField& p
) const
{
    const vector gHat =
        normalised
        (
            db().lookupObject<uniformDimensionedVectorField>("g").value()
        );

    const scalar h0 = - gHat & origin_;

    const vector UGasRefRel = UGasRef_ - UMean_->value(t);

    const scalar magUGasRefRel = mag(UGasRefRel);

    tmp<vectorField> tU = waveSuperposition::UGas(t, p);

    if (magUGasRefRel > 0)
    {
        atmBoundaryLayer atm
        (
            UGasRefRel/magUGasRefRel,
          - gHat,
            magUGasRefRel,
            h0 + hRef_,
            scalarField(p.size(), hWaveMax_ - hWaveMin_),
            scalarField(p.size(), h0 + hWaveMin_)
        );

        tU.ref() += atm.U(p);
    }

    return tU;
}


void Foam::waveAtmBoundaryLayerSuperposition::write(Ostream& os) const
{
    waveSuperposition::write(os);

    writeEntry(os, "UGasRef", UGasRef_);
    writeEntry(os, "hRef", hRef_);
    writeEntry(os, "hWaveMin", hWaveMin_);
    writeEntry(os, "hWaveMax", hWaveMax_);
}


// ************************************************************************* //
