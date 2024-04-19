/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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

#include "dynamicContactAngle.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace contactAngleModels
{
    defineTypeNameAndDebug(dynamic, 0);
    addToRunTimeSelectionTable(contactAngleModel, dynamic, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::contactAngleModels::dynamic::dynamic(const dictionary& dict)
:
    theta0_(dict.lookup<scalar>("theta0", unitDegrees)),
    uTheta_(dict.lookup<scalar>("uTheta", dimVelocity)),
    thetaAdv_(dict.lookup<scalar>("thetaAdv", unitDegrees)),
    thetaRec_(dict.lookup<scalar>("thetaRec", unitDegrees))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::contactAngleModels::dynamic::~dynamic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::contactAngleModels::dynamic::cosTheta
(
    const fvPatchVectorField& Up,
    const vectorField& nHat
) const
{
    const vectorField nf(Up.patch().nf());

    // Calculated the component of the velocity parallel to the wall
    vectorField U(Up.patchInternalField() - Up);
    U -= (nf & U)*nf;

    // Find the direction of the interface parallel to the wall and normalise
    vectorField n(nHat - (nf & nHat)*nf);
    n /= (mag(n) + small);

    // Calculate the velocity coefficient
    // from the component of U normal to the interface
    const scalarField uCoeff(tanh((n & U)/uTheta_));

    return cos
    (
        theta0_
      + thetaRec_ - theta0_*max(uCoeff, scalar(0))
      - thetaAdv_ - theta0_*min(uCoeff, scalar(0))
    );
}


void Foam::contactAngleModels::dynamic::write(Ostream& os) const
{
    writeEntry(os, "theta0", unitDegrees, theta0_);
    writeEntry(os, "uTheta", dimVelocity, uTheta_);
    writeEntry(os, "thetaAdv", unitDegrees, thetaAdv_);
    writeEntry(os, "thetaRec", unitDegrees, thetaRec_);
}


// ************************************************************************* //
