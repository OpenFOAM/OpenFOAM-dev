/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "phaseInterface.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::volScalarField& Foam::phaseInterface::alpha1() const
{
    return phase1();
}


inline const Foam::volScalarField& Foam::phaseInterface::alpha2() const
{
    return phase2();
}


inline const Foam::rhoFluidThermo& Foam::phaseInterface::thermo1() const
{
    return phase1().fluidThermo();
}


inline const Foam::rhoFluidThermo& Foam::phaseInterface::thermo2() const
{
    return phase2().fluidThermo();
}


inline const Foam::volScalarField& Foam::phaseInterface::rho1() const
{
    return phase1().rho();
}


inline const Foam::volScalarField& Foam::phaseInterface::rho2() const
{
    return phase2().rho();
}


inline const Foam::phaseModel& Foam::phaseInterface::phase1() const
{
    return phase1_;
}


inline const Foam::phaseModel& Foam::phaseInterface::phase2() const
{
    return phase2_;
}


inline bool Foam::phaseInterface::contains(const phaseModel& phase) const
{
    return &phase1_ == &phase || &phase2_ == &phase;
}


inline const Foam::phaseModel& Foam::phaseInterface::otherPhase
(
    const phaseModel& phase
) const
{
    if (&phase1_ == &phase)
    {
        return phase2_;
    }
    else if (&phase2_ == &phase)
    {
        return phase1_;
    }
    else
    {
        FatalErrorInFunction
            << "this phaseInterface does not contain phase " << phase.name()
            << exit(FatalError);

        return phase;
    }
}


inline Foam::label Foam::phaseInterface::index(const phaseModel& phase) const
{
    if (&phase1_ == &phase)
    {
        return 0;
    }
    else if (&phase2_ == &phase)
    {
        return 1;
    }
    else
    {
        FatalErrorInFunction
            << "this phaseInterface does not contain phase " << phase.name()
            << exit(FatalError);

        return -1;
    }
}


inline const Foam::phaseSystem& Foam::phaseInterface::fluid() const
{
    return phase1().fluid();
}


inline const Foam::fvMesh& Foam::phaseInterface::mesh() const
{
    return phase1().mesh();
}


inline const Foam::uniformDimensionedVectorField&
Foam::phaseInterface::g() const
{
    return g_;
}


// * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * * //

inline Foam::phaseInterface::const_iterator::const_iterator
(
    const phaseInterface& pair,
    const label index
)
:
    pair_(pair),
    index_(index)
{}


inline Foam::phaseInterface::const_iterator::const_iterator
(
    const phaseInterface& pair
)
:
    const_iterator(pair, 0)
{}


inline Foam::label Foam::phaseInterface::const_iterator::index() const
{
    return index_;
}


inline bool Foam::phaseInterface::const_iterator::operator==
(
    const const_iterator& iter
) const
{
    return (this->index_ == iter.index_);
}


inline bool Foam::phaseInterface::const_iterator::operator!=
(
    const const_iterator& iter
) const
{
    return !(this->operator==(iter));
}


inline const Foam::phaseModel&
Foam::phaseInterface::const_iterator::operator*() const
{
    if (index_ == 0)
    {
        return pair_.phase1_;
    }
    else
    {
        return pair_.phase2_;
    }
}


inline const Foam::phaseModel&
Foam::phaseInterface::const_iterator::operator()() const
{
    return operator*();
}


inline const Foam::phaseModel&
Foam::phaseInterface::const_iterator::otherPhase() const
{
    if (index_ == 0)
    {
        return pair_.phase2_;
    }
    else
    {
        return pair_.phase1_;
    }
}


inline Foam::phaseInterface::const_iterator&
Foam::phaseInterface::const_iterator::operator++()
{
    index_++;
    return *this;
}


inline Foam::phaseInterface::const_iterator
Foam::phaseInterface::const_iterator::operator++(int)
{
    const_iterator old = *this;
    this->operator++();
    return old;
}


inline Foam::phaseInterface::const_iterator Foam::phaseInterface::cbegin() const
{
    return const_iterator(*this);
}


inline Foam::phaseInterface::const_iterator Foam::phaseInterface::cend() const
{
    return const_iterator(*this, 2);
}


inline Foam::phaseInterface::const_iterator Foam::phaseInterface::begin() const
{
    return const_iterator(*this);
}


inline Foam::phaseInterface::const_iterator Foam::phaseInterface::end() const
{
    return const_iterator(*this, 2);
}


// ************************************************************************* //
