/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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

Class
    Foam::Identity

Description
    Templated identity and dual space identity tensors derived from
    SphericalTensor.

\*---------------------------------------------------------------------------*/

#ifndef Identity_H
#define Identity_H

#include "SphericalTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class Identity Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Identity
:
    public SphericalTensor<Cmpt>
{

public:

    //- Construct initialising the SphericalTensor to 1
    Identity()
    :
        SphericalTensor<Cmpt>(1)
    {}

    //- The identity type in the dual space
    class dual
    :
        public SphericalTensor<Cmpt>
    {
        //- Construct initialising the SphericalTensor to 1
        dual()
        :
            SphericalTensor<Cmpt>(1)
        {}
    };

    //- Return the identity in the dual space
    inline dual operator*()
    {
        return dual();
    }

    //- Return 1 for label
    inline explicit operator label() const
    {
        return 1;
    }

    //- Return 1 for scalar
    inline explicit operator scalar() const
    {
        return 1;
    }
};


// Global Identity tensor
static const Identity<scalar> I;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
