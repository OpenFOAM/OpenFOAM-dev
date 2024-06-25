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

#include "propellerDisk.H"
#include "fvMatrix.H"
#include "pimpleNoLoopControl.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::propellerDisk::addActuationDiskAxialInertialResistance
(
    vectorField& Usource,
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    const volVectorField::Internal& U
) const
{
    const labelUList& cells = set_.cells();
    const scalarField& V = mesh().V();

    volVectorField::Internal force
    (
        IOobject
        (
            typedName("force"),
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector(rho.dimensions()*U.dimensions()/dimTime, Zero)
    );

    const vector centre = diskCentre();
    const scalar delta = diskThickness(centre);

    // hub radius
    const scalar rHub = 0.5*dHub_;

    // propeller radius
    const scalar rProp = 0.5*dProp_;

    // Unit normal
    const vector nHat = normalised(diskNormal_);

    // Disk area
    const scalar A = (rProp - rHub)*(rProp - rHub)*pi;

    // Get reference to the n value
    const scalar& n = this->n();

    // Evaluate advance number J,
    // the volumetric mean advance speed of the disk
    const scalar J = this->J(U, nHat);

    // Calculate the mean velocity through the disk
    const scalar Udisk = n*dProp_*J;

    // ... and the volumetric flux through the disk
    const scalar phiDisk = Udisk*A;

    if (phiDisk > small)
    {
        // Correction of advance number J based on mometum theory
        // using the propeller curve Kt(J) and Kq(J)

        vector2D KtandKq(propellerFunction_->value(J));
        scalar Kt = KtandKq.x();
        scalar Kq = KtandKq.y();

        // Evaluate thrust based on Kt
        // Thrust units are force/rho = [m^4/s^2]
        scalar T = Kt*sqr(n)*pow4(dProp_);

        // Correction for speed based on momentum theory
        const scalar Ucorr = T/(2*phiDisk);
        const scalar Jcorr = (Udisk - Ucorr)/(n*dProp_);

        // Update Kt and Kq based on the corrected J value
        KtandKq = propellerFunction_->value(Jcorr);
        Kt = KtandKq.x();
        Kq = KtandKq.y();
        T = Kt*sqr(n)*pow4(dProp_);
        const scalar Q = Kq*sqr(n)*pow5(dProp_);

        // Correct the rotation speed from the current propulsion force
        correctn(T);

        // Distribute momentum source

        const scalar Ax =
            (105/8.0)*T/(delta*pi*(rProp - rHub)*(3*rHub + 4*rProp));
        const scalar At =
            (105/8.0)*Q/(delta*pi*rProp*(rProp - rHub)*(3*rHub + 4*rProp));

        forAll(cells, i)
        {
            const label celli = cells[i];

            const vector r = mesh().cellCentres()[celli] - centre;
            const scalar magr = mag(r);

            if (magr > rHub)
            {
                const scalar rPrime = magr/rProp;
                const scalar rHubPrime = rHub/rProp;

                // Normalised radius
                const scalar rStar =
                    min((rPrime - rHubPrime)/(1 - rHubPrime), 1);

                const scalar Faxial = Ax*rStar*sqrt(1 - rStar);
                const scalar Ftangential =
                    At*rStar*sqrt(1 - rStar)
                   /(rStar*(1 - rHubPrime) + rHubPrime);

                // A unit normal vector to the direction of the tangent
                const vector radiusOrtho = normalised(r ^ nHat);
                const vector rotationVector = radiusOrtho*rotationDir_;

                force[celli] =
                    alpha[celli]*rho[celli]
                   *(Faxial*nHat + Ftangential*rotationVector);

                Usource[celli] += V[celli]*force[celli];
            }
        }

        if
        (
            mesh().lookupObject<pimpleNoLoopControl>("solutionControl")
           .finalIter()
        )
        {
            if (logFile_.valid())
            {
                logFile_->writeTime(n, J, Jcorr, Udisk, Ucorr, Kt, Kq, T, Q);
            }

            if (mesh().time().writeTime())
            {
                force.write();
            }
        }
    }
}


// ************************************************************************* //
