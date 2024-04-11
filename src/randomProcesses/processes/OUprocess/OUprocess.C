/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "OUprocess.H"
#include "Kmesh.H"
#include "standardNormal.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

complexVector OUprocess::WeinerProcess(const scalar deltaT) const
{
    distributions::standardNormal stdNormal(rndGen_);

    return sqrt(deltaT)*complexVector
    (
        complex(stdNormal.sample(), stdNormal.sample()),
        complex(stdNormal.sample(), stdNormal.sample()),
        complex(stdNormal.sample(), stdNormal.sample())
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

OUprocess::OUprocess
(
    const Kmesh& kmesh,
    const scalar deltaT,
    const dictionary& OUdict
)
:
    rndGen_(label(0)),
    Kmesh_(kmesh),
    OUfield_(Kmesh_.size()),

    alpha_(OUdict.lookup<scalar>("alpha")),
    sigma_(OUdict.lookup<scalar>("sigma")),
    kUpper_(OUdict.lookup<scalar>("kUpper")),
    kLower_(OUdict.lookup<scalar>("kLower")),
    scale_((kUpper_ - kLower_)*pow(scalar(Kmesh_.size()), 1.0/vector::dim))
{
    const vectorField& K = Kmesh_;

    scalar sqrkUpper_ = sqr(kUpper_);
    scalar sqrkLower_ = sqr(kLower_) + small;
    scalar sqrK;

    forAll(OUfield_, i)
    {
        if ((sqrK = magSqr(K[i])) < sqrkUpper_ && sqrK > sqrkLower_)
        {
            OUfield_[i] = scale_*sigma_*WeinerProcess(deltaT);
        }
        else
        {
            OUfield_[i] = complexVector
            (
                complex(0, 0),
                complex(0, 0),
                complex(0, 0)
            );
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const complexVectorField& OUprocess::newField(const scalar deltaT) const
{
    const vectorField& K = Kmesh_;

    scalar sqrkUpper_ = sqr(kUpper_);
    scalar sqrkLower_ = sqr(kLower_) + small;
    scalar sqrK;

    forAll(OUfield_, i)
    {
        if ((sqrK = magSqr(K[i])) < sqrkUpper_ && sqrK > sqrkLower_)
        {
            OUfield_[i] =
                (1.0 - alpha_*deltaT)*OUfield_[i]
              + scale_*sigma_*WeinerProcess(deltaT);
        }
    }

    return OUfield_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
