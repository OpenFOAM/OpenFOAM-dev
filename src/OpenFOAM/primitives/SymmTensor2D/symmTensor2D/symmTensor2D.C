/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "symmTensor2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const symmTensor2D::typeName = "symmTensor2D";

template<>
const char* symmTensor2D::componentNames[] =
{
    "xx", "xy",
          "yy"
};

template<>
const symmTensor2D symmTensor2D::zero
(
    0, 0,
       0
);

template<>
const symmTensor2D symmTensor2D::one
(
    1, 1,
       1
);

template<>
const symmTensor2D symmTensor2D::max
(
    VGREAT, VGREAT,
            VGREAT
);

template<>
const symmTensor2D symmTensor2D::min
(
    -VGREAT, -VGREAT,
             -VGREAT
);

template<>
const symmTensor2D symmTensor2D::I
(
    1, 0,
       1
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
