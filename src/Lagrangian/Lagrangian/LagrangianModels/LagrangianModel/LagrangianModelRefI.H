/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LagrangianModelRef.H"
#include "error.H"
#include "nullObject.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::LagrangianModelRef::LagrangianModelRef()
:
    ref_(NullObjectRef<LagrangianModel>())
{}


inline Foam::LagrangianModelRef::LagrangianModelRef(const LagrangianModel& ref)
:
    ref_(ref)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline bool Foam::LagrangianModelRef::valid() const
{
    return notNull(ref_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline const Foam::LagrangianModel& Foam::LagrangianModelRef::operator()() const
{
    if (!valid())
    {
        FatalErrorInFunction
            << "Invalid Lagrangian model reference requested."
            << exit(FatalError);
    }

    return ref_;
}


// ************************************************************************* //
