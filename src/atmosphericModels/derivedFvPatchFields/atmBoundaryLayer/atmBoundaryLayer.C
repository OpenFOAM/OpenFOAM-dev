/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2024 OpenFOAM Foundation
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

#include "atmBoundaryLayer.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::atmBoundaryLayer::kappaDefault_ = 0.41;

const Foam::scalar Foam::atmBoundaryLayer::CmuDefault_ = 0.09;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::atmBoundaryLayer::init()
{
    if (mag(flowDir_) < small || mag(zDir_) < small)
    {
        FatalErrorInFunction
            << "magnitude of n or z must be greater than zero"
            << abort(FatalError);
    }

    // Ensure direction vectors are normalised
    flowDir_ /= mag(flowDir_);
    zDir_ /= mag(zDir_);

    Ustar_ = kappa_*Uref_/(log((Zref_ + z0_)/z0_));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::atmBoundaryLayer::atmBoundaryLayer()
:
    flowDir_(Zero),
    zDir_(Zero),
    kappa_(0.41),
    Cmu_(0.09),
    Uref_(0),
    Zref_(0),
    z0_(0),
    zGround_(0),
    Ustar_(0),
    offset_(false),
    Ulower_(0),
    kLower_(0),
    epsilonLower_(0)
{}


Foam::atmBoundaryLayer::atmBoundaryLayer
(
    const vector& flowDir,
    const vector& zDir,
    const scalar Uref,
    const scalar Zref,
    const scalarField& z0,
    const scalarField& zGround,
    const scalar kappa,
    const scalar Cmu,
    const scalar Ulower,
    const scalar kLower,
    const scalar epsilonLower
)
:
    flowDir_(flowDir),
    zDir_(zDir),
    kappa_(kappa),
    Cmu_(Cmu),
    Uref_(Uref),
    Zref_(Zref),
    z0_(z0),
    zGround_(zGround),
    Ustar_(z0.size()),
    offset_(Ulower != 0),
    Ulower_(Ulower),
    kLower_(kLower),
    epsilonLower_(epsilonLower)
{
    init();
}


Foam::atmBoundaryLayer::atmBoundaryLayer
(
    const vectorField& p,
    const dictionary& dict
)
:
    flowDir_(dict.lookup<vector>("flowDir", dimless)),
    zDir_(dict.lookup<vector>("zDir", dimless)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", dimless, kappaDefault_)),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", dimless, CmuDefault_)),
    Uref_(dict.lookup<scalar>("Uref", dimVelocity)),
    Zref_(dict.lookup<scalar>("Zref", dimLength)),
    z0_("z0", dimLength, dict, p.size()),
    zGround_("zGround", dimLength, dict, p.size()),
    Ustar_(p.size()),
    offset_(dict.found("Ulower")),
    Ulower_(dict.lookupOrDefault<scalar>("Ulower", dimVelocity, 0)),
    kLower_(dict.lookupOrDefault<scalar>("kLower", dimEnergy/dimMass, 0)),
    epsilonLower_
    (
        dict.lookupOrDefault<scalar>
        (
            "epsilonLower",
            dimEnergy/dimMass/dimTime,
            0
        )
    )
{
    init();
}


Foam::atmBoundaryLayer::atmBoundaryLayer
(
    const atmBoundaryLayer& abl,
    const fieldMapper& mapper
)
:
    flowDir_(abl.flowDir_),
    zDir_(abl.zDir_),
    kappa_(abl.kappa_),
    Cmu_(abl.Cmu_),
    Uref_(abl.Uref_),
    Zref_(abl.Zref_),
    z0_(mapper(abl.z0_)),
    zGround_(mapper(abl.zGround_)),
    Ustar_(mapper(abl.Ustar_)),
    offset_(abl.offset_),
    Ulower_(abl.Ulower_),
    kLower_(abl.kLower_),
    epsilonLower_(abl.epsilonLower_)
{}


Foam::atmBoundaryLayer::atmBoundaryLayer(const atmBoundaryLayer& abl)
:
    flowDir_(abl.flowDir_),
    zDir_(abl.zDir_),
    kappa_(abl.kappa_),
    Cmu_(abl.Cmu_),
    Uref_(abl.Uref_),
    Zref_(abl.Zref_),
    z0_(abl.z0_),
    zGround_(abl.zGround_),
    Ustar_(abl.Ustar_),
    offset_(abl.offset_),
    Ulower_(abl.Ulower_),
    kLower_(abl.kLower_),
    epsilonLower_(abl.epsilonLower_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::atmBoundaryLayer::map
(
    const atmBoundaryLayer& blptf,
    const fieldMapper& mapper
)
{
    mapper(z0_, blptf.z0_);
    mapper(zGround_, blptf.zGround_);
    mapper(Ustar_, blptf.Ustar_);
}


void Foam::atmBoundaryLayer::reset(const atmBoundaryLayer& blptf)
{
    z0_.reset(blptf.z0_);
    zGround_.reset(blptf.zGround_);
    Ustar_.reset(blptf.Ustar_);
}


Foam::tmp<Foam::vectorField> Foam::atmBoundaryLayer::U
(
    const vectorField& p
) const
{
    const scalarField Un
    (
        (Ustar_/kappa_)
       *log(max((zDir_ & p) - zGround_ + z0_, z0_)/z0_)
    );

    if (offset_)
    {
        return flowDir_*Un + flowDir_*Ulower_;
    }
    else
    {
        return flowDir_*Un;
    }
}


Foam::tmp<Foam::scalarField> Foam::atmBoundaryLayer::k
(
    const vectorField& p
) const
{
    tmp<scalarField> tk
    (
        sqr(Ustar_)/sqrt(Cmu_)
    );

    if (offset_)
    {
        const scalarField z((zDir_ & p) - zGround_);
        tk.ref() = pos0(z)*tk() + neg(z)*kLower_;
    }

    return tk;
}


Foam::tmp<Foam::scalarField> Foam::atmBoundaryLayer::epsilon
(
    const vectorField& p
) const
{
    tmp<scalarField> tepsilon
    (
        pow3(Ustar_)/(kappa_*((zDir_ & p) - zGround_ + z0_))
    );

    if (offset_)
    {
        const scalarField z((zDir_ & p) - zGround_);
        tepsilon.ref() = pos0(z)*tepsilon() + neg(z)*epsilonLower_;
    }

    return tepsilon;
}


void Foam::atmBoundaryLayer::write(Ostream& os) const
{
    writeEntry(os, "z0", z0_) ;
    writeEntry(os, "flowDir", flowDir_);
    writeEntry(os, "zDir", zDir_);
    writeEntry(os, "kappa", kappa_);
    writeEntry(os, "Cmu", Cmu_);
    writeEntry(os, "Uref", Uref_);
    writeEntry(os, "Zref", Zref_);

    if (offset_)
    {
        writeEntry(os, "Ulower", Ulower_);
        writeEntry(os, "kLower", kLower_);
        writeEntry(os, "epsilonLower", epsilonLower_);
    }

    writeEntry(os, "zGround", zGround_);
}


// ************************************************************************* //
