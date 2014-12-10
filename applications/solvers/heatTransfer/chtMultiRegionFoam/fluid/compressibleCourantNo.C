/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "compressibleCourantNo.H"
#include "fvc.H"

Foam::scalar Foam::compressibleCourantNo
(
    const fvMesh& mesh,
    const Time& runTime,
    const volScalarField& rho,
    const surfaceScalarField& phi
)
{
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().internalField()
      / rho.internalField()
    );

    scalar CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    scalar meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();

    Info<< "Region: " << mesh.name() << " Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;

    return CoNum;
}


// ************************************************************************* //
