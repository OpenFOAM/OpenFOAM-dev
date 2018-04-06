/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "moleculeCloud.H"
#include "molecule.H"
#include "Random.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tensor Foam::molecule::rotationTensorX(scalar phi) const
{
    return tensor
    (
        1, 0, 0,
        0, Foam::cos(phi), -Foam::sin(phi),
        0, Foam::sin(phi), Foam::cos(phi)
    );
}


Foam::tensor Foam::molecule::rotationTensorY(scalar phi) const
{
    return tensor
    (
        Foam::cos(phi), 0, Foam::sin(phi),
        0, 1, 0,
        -Foam::sin(phi), 0, Foam::cos(phi)
    );
}


Foam::tensor Foam::molecule::rotationTensorZ(scalar phi) const
{
    return tensor
    (
        Foam::cos(phi), -Foam::sin(phi), 0,
        Foam::sin(phi), Foam::cos(phi), 0,
        0, 0, 1
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::molecule::move
(
    moleculeCloud& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    const constantProperties& constProps(cloud.constProps(id_));

    if (td.part() == 0)
    {
        // First leapfrog velocity adjust part, required before tracking+force
        // part

        v_ += 0.5*trackTime*a_;

        pi_ += 0.5*trackTime*tau_;
    }
    else if (td.part() == 1)
    {
        // Leapfrog tracking part

        while (td.keepParticle && !td.switchProcessor && stepFraction() < 1)
        {
            const scalar f = 1 - stepFraction();
            trackToAndHitFace(f*trackTime*v_, f, cloud, td);
        }
    }
    else if (td.part() == 2)
    {
        // Leapfrog orientation adjustment, carried out before force calculation
        // but after tracking stage, i.e. rotation carried once linear motion
        // complete.

        if (!constProps.pointMolecule())
        {
            const diagTensor& momentOfInertia(constProps.momentOfInertia());

            tensor R;

            if (!constProps.linearMolecule())
            {
                R = rotationTensorX(0.5*trackTime*pi_.x()/momentOfInertia.xx());
                pi_ = pi_ & R;
                Q_ = Q_ & R;
            }

            R = rotationTensorY(0.5*trackTime*pi_.y()/momentOfInertia.yy());
            pi_ = pi_ & R;
            Q_ = Q_ & R;

            R = rotationTensorZ(trackTime*pi_.z()/momentOfInertia.zz());
            pi_ = pi_ & R;
            Q_ = Q_ & R;

            R = rotationTensorY(0.5*trackTime*pi_.y()/momentOfInertia.yy());
            pi_ = pi_ & R;
            Q_ = Q_ & R;

            if (!constProps.linearMolecule())
            {
                R = rotationTensorX(0.5*trackTime*pi_.x()/momentOfInertia.xx());
                pi_ = pi_ & R;
                Q_ = Q_ & R;
            }
        }

        setSitePositions(constProps);
    }
    else if (td.part() == 3)
    {
        // Second leapfrog velocity adjust part, required after tracking+force
        // part

        scalar m = constProps.mass();

        a_ = Zero;

        tau_ = Zero;

        forAll(siteForces_, s)
        {
            const vector& f = siteForces_[s];

            a_ += f/m;

            tau_ += (constProps.siteReferencePositions()[s] ^ (Q_.T() & f));
        }

        v_ += 0.5*trackTime*a_;

        pi_ += 0.5*trackTime*tau_;

        if (constProps.pointMolecule())
        {
            tau_ = Zero;

            pi_ = Zero;
        }

        if (constProps.linearMolecule())
        {
            tau_.x() = 0.0;

            pi_.x() = 0.0;
        }
    }
    else
    {
        FatalErrorInFunction
            << td.part() << " is an invalid part of the integration method."
            << abort(FatalError);
    }

    return td.keepParticle;
}


void Foam::molecule::transformProperties(const tensor& T)
{
    particle::transformProperties(T);

    Q_ = T & Q_;

    v_ = transform(T, v_);

    a_ = transform(T, a_);

    pi_ = Q_.T() & transform(T, Q_ & pi_);

    tau_ = Q_.T() & transform(T, Q_ & tau_);

    rf_ = transform(T, rf_);

    sitePositions_ = position() + (T & (sitePositions_ - position()));

    siteForces_ = T & siteForces_;
}


void Foam::molecule::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);

    if (special_ == SPECIAL_TETHERED)
    {
        specialPosition_ += separation;
    }

    sitePositions_ = sitePositions_ + separation;
}


void Foam::molecule::setSitePositions(const constantProperties& constProps)
{
    sitePositions_ = position() + (Q_ & constProps.siteReferencePositions());
}


void Foam::molecule::setSiteSizes(label size)
{
    sitePositions_.setSize(size);

    siteForces_.setSize(size);
}


bool Foam::molecule::hitPatch(moleculeCloud&, trackingData&)
{
    return false;
}


void Foam::molecule::hitProcessorPatch(moleculeCloud&, trackingData& td)
{
    td.switchProcessor = true;
}


void Foam::molecule::hitWallPatch(moleculeCloud&, trackingData&)
{
    const vector nw = normal();
    const scalar vn = v_ & nw;

    // Specular reflection
    if (vn > 0)
    {
        v_ -= 2*vn*nw;
    }
}


// ************************************************************************* //
