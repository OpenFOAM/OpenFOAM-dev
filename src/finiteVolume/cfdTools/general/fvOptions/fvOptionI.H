/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::fv::option::name() const
{
    return name_;
}


inline const Foam::fvMesh& Foam::fv::option::mesh() const
{
    return mesh_;
}


inline const Foam::dictionary& Foam::fv::option::coeffs() const
{
    return coeffs_;
}


inline bool Foam::fv::option::active() const
{
    return active_;
}


inline void Foam::fv::option::setApplied(const label fieldi)
{
    applied_[fieldi] = true;
}


inline Foam::Switch& Foam::fv::option::active()
{
    return active_;
}


// ************************************************************************* //
