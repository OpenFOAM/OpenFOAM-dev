/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::LagrangianModelRef::LagrangianModelRef()
:
    source_(NullObjectRef<LagrangianSource>()),
    S_(NullObjectRef<LagrangianSubScalarField>()),
    injection_(NullObjectRef<LagrangianInjection>())
{}


inline Foam::LagrangianModelRef::LagrangianModelRef
(
    const LagrangianSource& source,
    const LagrangianSubScalarField& S
)
:
    source_(source),
    S_(S),
    injection_(NullObjectRef<LagrangianInjection>())
{}


inline Foam::LagrangianModelRef::LagrangianModelRef
(
    const LagrangianInjection& injection
)
:
    source_(NullObjectRef<LagrangianSource>()),
    S_(NullObjectRef<LagrangianSubScalarField>()),
    injection_(injection)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline bool Foam::LagrangianModelRef::isSource() const
{
    return notNull(source_);
}


inline const Foam::LagrangianSource&
Foam::LagrangianModelRef::source() const
{
    if (!isSource())
    {
        FatalErrorInFunction
            << "Invalid Lagrangian source reference requested."
            << exit(FatalError);
    }

    return source_;
}


inline const Foam::LagrangianSubScalarField&
Foam::LagrangianModelRef::S() const
{
    if (!isSource())
    {
        FatalErrorInFunction
            << "Invalid Lagrangian source reference requested."
            << exit(FatalError);
    }

    return S_;
}


inline bool Foam::LagrangianModelRef::isInjection() const
{
    return notNull(injection_);
}


inline const Foam::LagrangianInjection&
Foam::LagrangianModelRef::injection() const
{
    if (!isInjection())
    {
        FatalErrorInFunction
            << "Invalid Lagrangian injection reference requested."
            << exit(FatalError);
    }

    return injection_;
}


// ************************************************************************* //
