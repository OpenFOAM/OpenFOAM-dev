/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "correctContactAngle.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::correctContactAngle
(
    const volScalarField& a1,
    const volScalarField& a2,
    const volVectorField::Boundary& Ubf,
    const dimensionedScalar& deltaN,
    surfaceVectorField::Boundary& nHatbf
)
{
    const volScalarField::Boundary& a1bf = a1.boundaryField();
    const volScalarField::Boundary& a2bf = a2.boundaryField();

    const fvBoundaryMesh& boundary = a1.mesh().boundary();

    forAll(boundary, patchi)
    {
        const fvPatchScalarField& a1p = a1bf[patchi];
        const fvPatchScalarField& a2p = a2bf[patchi];

        const bool a1pIsCa = isA<alphaContactAngleFvPatchScalarField>(a1p);
        const bool a2pIsCa = isA<alphaContactAngleFvPatchScalarField>(a2p);

        if (a1pIsCa || a2pIsCa)
        {
            if (a1pIsCa != a2pIsCa)
            {
                FatalErrorInFunction
                    << alphaContactAngleFvPatchScalarField::typeName
                    << " boundary condition specified on patch "
                    << boundary[patchi].name() << " for "
                    << (a1pIsCa ? a1.name() : a2.name()) << " but not for "
                    << (a2pIsCa ? a1.name() : a2.name())
                    << exit(FatalError);
            }

            const alphaContactAngleFvPatchScalarField& a1ca =
                refCast<const alphaContactAngleFvPatchScalarField>(a1p);
            const alphaContactAngleFvPatchScalarField& a2ca =
                refCast<const alphaContactAngleFvPatchScalarField>(a2p);

            const bool a1caHasProps = a1ca.thetaProps().found(a2.group());
            const bool a2caHasProps = a2ca.thetaProps().found(a1.group());

            if (!a1caHasProps && !a2caHasProps)
            {
                FatalErrorInFunction
                    << "Neither "
                    << alphaContactAngleFvPatchScalarField::typeName
                    << " boundary condition specified on patch "
                    << boundary[patchi].name()
                    << " for " << a1.name() << " and " << a2.name()
                    << " contains properties for the other phase"
                    << exit(FatalError);
            }

            if
            (
                a1caHasProps && a2caHasProps
             && a1ca.thetaProps()[a2.group()]
             != a2ca.thetaProps()[a1.group()].reversed()
            )
            {
                FatalErrorInFunction
                    << "The "
                    << alphaContactAngleFvPatchScalarField::typeName
                    << " boundary conditions specified on patch "
                    << boundary[patchi].name()
                    << " for " << a1.name() << " and " << a2.name()
                    << " contain inconsistent properties"
                    << exit(FatalError);
            }

            const alphaContactAngleFvPatchScalarField::contactAngleProperties
                tp = a1caHasProps
              ? a1ca.thetaProps()[a2.group()]
              : a2ca.thetaProps()[a1.group()].reversed();

            const vectorField np(a1.mesh().boundary()[patchi].nf());

            vectorField& nHatp = nHatbf[patchi];

            // Calculate the contact angle
            scalarField theta(np.size(), degToRad(tp.theta0()));

            // Calculate the dynamic contact angle if required
            if (tp.dynamic())
            {
                const scalar uTheta = tp.uTheta();
                const scalar thetaA = degToRad(tp.thetaA());
                const scalar thetaR = degToRad(tp.thetaR());

                // Calculated the component of the velocity parallel to the wall
                vectorField Uwall
                (
                    Ubf[patchi].patchInternalField() - Ubf[patchi]
                );
                Uwall -= (np & Uwall)*np;

                // Find the direction of the interface parallel to the wall
                vectorField nWall(nHatp - (np & nHatp)*np);
                nWall /= (mag(nWall) + small);

                // Calculate Uwall resolved normal to the interface parallel to
                // the interface
                const scalarField uwall(nWall & Uwall);

                theta += (thetaA - thetaR)*tanh(uwall/uTheta);
            }

            // Reset nHatp to correspond to the contact angle
            const scalarField a12(nHatp & np);
            const scalarField b1(cos(theta));
            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }
            const scalarField det(1 - a12*a12);
            const scalarField a((b1 - a12*b2)/det);
            const scalarField b((b2 - a12*b1)/det);

            nHatp = a*np + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN.value());
        }
    }
}


// ************************************************************************* //
