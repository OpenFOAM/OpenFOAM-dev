/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "skewCorrectionVectors.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(skewCorrectionVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::skewCorrectionVectors::skewCorrectionVectors(const fvMesh& mesh)
:
    DemandDrivenMeshObject
    <
        fvMesh,
        MoveableMeshObject,
        skewCorrectionVectors
    >(mesh),
    skew_(false),
    skewCorrectionVectors_
    (
        IOobject
        (
            "skewCorrectionVectors",
            mesh.pointsInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimless
    )
{
    calcSkewCorrectionVectors();
}


Foam::skewCorrectionVectors::~skewCorrectionVectors()
{}


void Foam::skewCorrectionVectors::calcSkewCorrectionVectors()
{
    if (debug)
    {
        InfoInFunction << "Calculating skew correction vectors" << endl;
    }

    // Set local references to mesh data
    const volVectorField& C = mesh().C();
    const surfaceVectorField& Cf = mesh().Cf();
    const surfaceVectorField& Sf = mesh().Sf();

    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        vector d = C[nei] - C[own];
        vector Cpf = Cf[facei] - C[own];

        skewCorrectionVectors_[facei] =
            Cpf - ((Sf[facei] & Cpf)/(Sf[facei] & d))*d;
    }

    surfaceVectorField::Boundary& skewCorrVecsBf =
        skewCorrectionVectors_.boundaryFieldRef();

    forAll(skewCorrVecsBf, patchi)
    {
        fvsPatchVectorField& patchSkewCorrVecs = skewCorrVecsBf[patchi];

        if (!patchSkewCorrVecs.coupled())
        {
            patchSkewCorrVecs = Zero;
        }
        else
        {
            const fvPatch& p = patchSkewCorrVecs.patch();
            const labelUList& faceCells = p.faceCells();
            const vectorField& patchFaceCentres = Cf.boundaryField()[patchi];
            const vectorField& patchSf = Sf.boundaryField()[patchi];
            const vectorField patchD(p.delta());

            forAll(p, patchFacei)
            {
                vector Cpf =
                    patchFaceCentres[patchFacei] - C[faceCells[patchFacei]];

                patchSkewCorrVecs[patchFacei] =
                    Cpf
                  - (
                        (patchSf[patchFacei] & Cpf)/
                        (patchSf[patchFacei] & patchD[patchFacei])
                    )*patchD[patchFacei];
            }
        }
    }

    scalar skewCoeff = 0.0;

    if (Sf.primitiveField().size())
    {
        skewCoeff =
            max(mag(skewCorrectionVectors_)*mesh().deltaCoeffs()).value();
    }

    if (debug)
    {
        InfoInFunction << "skew coefficient = " << skewCoeff << endl;
    }

    if (skewCoeff < 1e-5)
    {
        skew_ = false;
    }
    else
    {
        skew_ = true;
    }

    if (debug)
    {
        Info<< "    Finished constructing skew correction vectors" << endl;
    }
}


bool Foam::skewCorrectionVectors::movePoints()
{
    calcSkewCorrectionVectors();
    return true;
}


// ************************************************************************* //
