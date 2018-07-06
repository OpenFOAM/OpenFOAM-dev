/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::CompactSpatialTensorT<Cmpt>::CompactSpatialTensorT()
{}


template<class Cmpt>
inline Foam::CompactSpatialTensorT<Cmpt>::CompactSpatialTensorT
(
    const Foam::zero
)
:
    CompactSpatialTensorT::msType(Zero)
{}


template<class Cmpt>
inline Foam::CompactSpatialTensorT<Cmpt>::CompactSpatialTensorT
(
    const typename CompactSpatialTensorT::msType& ms
)
:
    CompactSpatialTensorT::msType(ms)
{}


template<class Cmpt>
inline Foam::CompactSpatialTensorT<Cmpt>::CompactSpatialTensorT
(
    const Cmpt& t00, const Cmpt& t01, const Cmpt& t02,
    const Cmpt& t10, const Cmpt& t11, const Cmpt& t12,
    const Cmpt& t20, const Cmpt& t21, const Cmpt& t22,
    const Cmpt& t30, const Cmpt& t31, const Cmpt& t32,
    const Cmpt& t40, const Cmpt& t41, const Cmpt& t42,
    const Cmpt& t50, const Cmpt& t51, const Cmpt& t52
)
{
    this->v_[0] = t00;
    this->v_[1] = t01;
    this->v_[2] = t02;

    this->v_[3 + 0] = t10;
    this->v_[3 + 1] = t11;
    this->v_[3 + 2] = t12;

    this->v_[6 + 0] = t20;
    this->v_[6 + 1] = t21;
    this->v_[6 + 2] = t22;

    this->v_[9 + 0] = t30;
    this->v_[9 + 1] = t31;
    this->v_[9 + 2] = t32;

    this->v_[12 + 0] = t40;
    this->v_[12 + 1] = t41;
    this->v_[12 + 2] = t42;

    this->v_[15 + 0] = t50;
    this->v_[15 + 1] = t51;
    this->v_[15 + 2] = t52;
}


template<class Cmpt>
inline Foam::CompactSpatialTensorT<Cmpt>::CompactSpatialTensorT(Istream& is)
:
    CompactSpatialTensorT::msType(is)
{}


// ************************************************************************* //
