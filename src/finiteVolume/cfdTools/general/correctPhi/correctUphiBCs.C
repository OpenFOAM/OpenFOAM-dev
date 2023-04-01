/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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

#include "fvCorrectPhi.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::correctUphiBCs
(
    volVectorField& U,
    surfaceScalarField& phi,
    const bool evaluateUBCs
)
{
    const fvMesh& mesh = U.mesh();

    if (mesh.changing())
    {
        volVectorField::Boundary& Ubf = U.boundaryFieldRef();
        surfaceScalarField::Boundary& phibf = phi.boundaryFieldRef();

        if (evaluateUBCs)
        {
            forAll(Ubf, patchi)
            {
                if (Ubf[patchi].fixesValue())
                {
                    Ubf[patchi].initEvaluate();
                }
            }
        }

        forAll(Ubf, patchi)
        {
            if (Ubf[patchi].fixesValue())
            {
                if (evaluateUBCs)
                {
                    Ubf[patchi].evaluate();
                }

                phibf[patchi] = Ubf[patchi] & mesh.Sf().boundaryField()[patchi];
            }
        }
    }
}


void Foam::correctUphiBCs
(
    const volScalarField& rho,
    volVectorField& U,
    surfaceScalarField& phi,
    const bool evaluateUBCs
)
{
    const fvMesh& mesh = U.mesh();

    if (mesh.changing())
    {
        volVectorField::Boundary& Ubf = U.boundaryFieldRef();
        surfaceScalarField::Boundary& phibf = phi.boundaryFieldRef();

        if (evaluateUBCs)
        {
            forAll(Ubf, patchi)
            {
                if (Ubf[patchi].fixesValue())
                {
                    Ubf[patchi].initEvaluate();
                }
            }
        }

        forAll(Ubf, patchi)
        {
            if (Ubf[patchi].fixesValue())
            {
                if (evaluateUBCs)
                {
                    Ubf[patchi].evaluate();
                }

                phibf[patchi] =
                    rho.boundaryField()[patchi]
                   *(
                        Ubf[patchi]
                      & mesh.Sf().boundaryField()[patchi]
                    );
            }
        }
    }
}


// ************************************************************************* //
