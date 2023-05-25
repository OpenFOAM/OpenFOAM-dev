/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "isothermalFilm.H"
#include "fvcGrad.H"
#include "filmContactAngleFvPatchScalarField.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField::Internal>
Foam::solvers::isothermalFilm::contactForce(const volScalarField& sigma) const
{
    tmp<volVectorField::Internal> tForce
    (
        volVectorField::Internal::New
        (
            typedName("contactForce"),
            mesh,
            dimensionedVector(dimForce/dimVolume, Zero)
        )
    );

    vectorField& force = tForce.ref();

    const pointField& points = mesh.points();
    const labelUList& own = mesh.owner();
    const labelUList& nbr = mesh.neighbour();

    const scalarField& V = mesh.V();

    const volVectorField gradDelta(fvc::grad(constrainedField(delta)));

    const scalar deltaWet = this->deltaWet.value();

    boolList contactCells(delta.size(), false);

    forAll(nbr, facei)
    {
        const label cello = own[facei];
        const label celln = nbr[facei];

        if ((delta[cello] > deltaWet) && (delta[celln] < deltaWet))
        {
            contactCells[cello] = true;
        }
        else if ((delta[cello] < deltaWet) && (delta[celln] > deltaWet))
        {
            contactCells[celln] = true;
        }
    }

    // Filter for film wall and surface patches
    labelHashSet wallAndSurfacePatches(wallPatchIDs);
    wallAndSurfacePatches.insert(surfacePatchID);

    const volScalarField::Boundary& deltaBf = delta.boundaryField();

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p(mesh.boundary()[patchi]);

        // For internal coupled patches
        if (p.coupled() && !wallAndSurfacePatches.found(patchi))
        {
            tmp<scalarField> tdeltan = deltaBf[patchi].patchNeighbourField();
            const scalarField& deltan = tdeltan();

            forAll(deltan, facei)
            {
                const label celli = p.faceCells()[facei];

                if ((delta[celli] > deltaWet) && (deltan[facei] < deltaWet))
                {
                    contactCells[celli] = true;
                }
            }
        }
    }

    forAll(wallPatchIDs, i)
    {
        const label patchi = wallPatchIDs[i];

        if
        (
            isA<filmContactAngleFvPatchScalarField>
            (
                delta.boundaryField()[patchi]
            )
        )
        {
            const filmContactAngleFvPatchScalarField& contactAngle
            (
                refCast<const filmContactAngleFvPatchScalarField>
                (
                    delta.boundaryField()[patchi]
                )
            );

            const vectorField np
            (
                gradDelta.boundaryField()[patchi]
               /(mag(gradDelta.boundaryField()[patchi]) + rootVSmall)
            );

            const scalarField cosTheta
            (
                contactAngle.cosTheta(U.boundaryField()[patchi], np)
            );

            const fvPatch& p(mesh.boundary()[patchi]);

            forAll(p, facei)
            {
                const label celli = p.faceCells()[facei];

                if (contactCells[celli])
                {
                    const vector& n = np[facei];

                    // Estimate the length of the contact line across the face
                    scalar clen = 0;
                    const face& f(p.patch()[facei]);
                    forAll(f, i)
                    {
                        clen += mag
                        (
                            n ^ (points[f[f.fcIndex(i)]] - points[f[i]])
                        );
                    }
                    clen /= 2;

                    // Calculate the contact force for the film wall face
                    force[celli] =
                        n*sigma[celli]*(1 - cosTheta[facei])*clen/V[celli];
                }
            }
        }
    }

    return tForce;
}


// ************************************************************************* //
