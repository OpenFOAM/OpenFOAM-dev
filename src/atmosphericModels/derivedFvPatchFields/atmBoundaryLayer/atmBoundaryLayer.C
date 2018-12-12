/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmBoundaryLayer::atmBoundaryLayer()
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


atmBoundaryLayer::atmBoundaryLayer(const vectorField& p, const dictionary& dict)
:
    flowDir_(dict.lookup("flowDir")),
    zDir_(dict.lookup("zDir")),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    Uref_(readScalar(dict.lookup("Uref"))),
    Zref_(readScalar(dict.lookup("Zref"))),
    z0_("z0", dict, p.size()),
    zGround_("zGround", dict, p.size()),
    Ustar_(p.size()),
    offset_(dict.found("Ulower")),
    Ulower_(dict.lookupOrDefault<scalar>("Ulower", 0)),
    kLower_(dict.lookupOrDefault<scalar>("kLower", 0)),
    epsilonLower_(dict.lookupOrDefault<scalar>("epsilonLower", 0))
{
    if (mag(flowDir_) < small || mag(zDir_) < small)
    {
        FatalErrorInFunction
            << "magnitude of n or z must be greater than zero"
            << abort(FatalError);
    }

    // Ensure direction vectors are normalized
    flowDir_ /= mag(flowDir_);
    zDir_ /= mag(zDir_);

    Ustar_ = kappa_*Uref_/(log((Zref_ + z0_)/z0_));
}


atmBoundaryLayer::atmBoundaryLayer
(
    const atmBoundaryLayer& abl,
    const fvPatchFieldMapper& mapper
)
:
    flowDir_(abl.flowDir_),
    zDir_(abl.zDir_),
    kappa_(abl.kappa_),
    Cmu_(abl.Cmu_),
    Uref_(abl.Uref_),
    Zref_(abl.Zref_),
    z0_(abl.z0_, mapper),
    zGround_(abl.zGround_, mapper),
    Ustar_(abl.Ustar_, mapper),
    offset_(abl.offset_),
    Ulower_(abl.Ulower_),
    kLower_(abl.kLower_),
    epsilonLower_(abl.epsilonLower_)
{}


atmBoundaryLayer::atmBoundaryLayer(const atmBoundaryLayer& abl)
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

void atmBoundaryLayer::autoMap(const fvPatchFieldMapper& m)
{
    z0_.autoMap(m);
    zGround_.autoMap(m);
    Ustar_.autoMap(m);
}


void atmBoundaryLayer::rmap
(
    const atmBoundaryLayer& blptf,
    const labelList& addr
)
{
    z0_.rmap(blptf.z0_, addr);
    zGround_.rmap(blptf.zGround_, addr);
    Ustar_.rmap(blptf.Ustar_, addr);
}


tmp<vectorField> atmBoundaryLayer::U(const vectorField& p) const
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


tmp<scalarField> atmBoundaryLayer::k(const vectorField& p) const
{
    tmp<scalarField> tk
    (
        sqr(Ustar_)/sqrt(Cmu_)
    );

    if (offset_)
    {
        const scalarField z = (zDir_ & p) - zGround_;
        tk.ref() = pos0(z)*tk() + neg(z)*kLower_;
    }

    return tk;
}


tmp<scalarField> atmBoundaryLayer::epsilon(const vectorField& p) const
{
    tmp<scalarField> tepsilon
    (
        pow3(Ustar_)/(kappa_*((zDir_ & p) - zGround_ + z0_))
    );

    if (offset_)
    {
        const scalarField z = (zDir_ & p) - zGround_;
        tepsilon.ref() = pos0(z)*tepsilon() + neg(z)*epsilonLower_;
    }

    return tepsilon;
}


void atmBoundaryLayer::write(Ostream& os) const
{
    z0_.writeEntry("z0", os) ;
    os.writeKeyword("flowDir")
        << flowDir_ << token::END_STATEMENT << nl;
    os.writeKeyword("zDir")
        << zDir_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa")
        << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cmu")
        << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("Uref")
        << Uref_ << token::END_STATEMENT << nl;
    os.writeKeyword("Zref")
        << Zref_ << token::END_STATEMENT << nl;

    if (offset_)
    {
        os.writeKeyword("Ulower")
            << Ulower_ << token::END_STATEMENT << nl;
        os.writeKeyword("kLower")
            << kLower_ << token::END_STATEMENT << nl;
        os.writeKeyword("epsilonLower")
            << epsilonLower_ << token::END_STATEMENT << nl;
    }

    zGround_.writeEntry("zGround", os) ;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
