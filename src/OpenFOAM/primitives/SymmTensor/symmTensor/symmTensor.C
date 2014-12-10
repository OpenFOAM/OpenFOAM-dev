/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "symmTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const symmTensor::typeName = "symmTensor";

template<>
const char* symmTensor::componentNames[] =
{
    "xx", "xy", "xz",
          "yy", "yz",
                "zz"
};

template<>
const symmTensor symmTensor::zero
(
    0, 0, 0,
       0, 0,
          0
);

template<>
const symmTensor symmTensor::one
(
    1, 1, 1,
       1, 1,
          1
);

template<>
const symmTensor symmTensor::max
(
    VGREAT, VGREAT, VGREAT,
            VGREAT, VGREAT,
                    VGREAT
);

template<>
const symmTensor symmTensor::min
(
    -VGREAT, -VGREAT, -VGREAT,
             -VGREAT, -VGREAT,
                      -VGREAT
);

template<>
const symmTensor symmTensor::I
(
    1, 0, 0,
       1, 0,
          1
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
