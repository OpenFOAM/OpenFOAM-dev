/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2020 OpenFOAM Foundation
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

inline const Foam::phaseModel& Foam::phasePair::phase1() const
{
    return phase1_;
}


inline const Foam::phaseModel& Foam::phasePair::phase2() const
{
    return phase2_;
}


inline bool Foam::phasePair::contains(const phaseModel& phase) const
{
    return &phase1_ == &phase || & phase2_ == &phase;
}


inline const Foam::phaseModel& Foam::phasePair::otherPhase
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
            << "this phasePair does not contain phase " << phase.name()
            << exit(FatalError);

        return phase;
    }
}


inline Foam::label Foam::phasePair::index(const phaseModel& phase) const
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
            << "this phasePair does not contain phase " << phase.name()
            << exit(FatalError);

        return -1;
    }
}


inline const Foam::uniformDimensionedVectorField& Foam::phasePair::g() const
{
    return g_;
}


inline const Foam::fvMesh& Foam::phasePair::mesh() const
{
    return phase1_.mesh();
}


// * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * * //

inline Foam::phasePair::const_iterator::const_iterator
(
    const phasePair& pair,
    const label index
)
:
    pair_(pair),
    index_(index)
{}


inline Foam::phasePair::const_iterator::const_iterator(const phasePair& pair)
:
    const_iterator(pair, 0)
{}


inline Foam::label Foam::phasePair::const_iterator::index() const
{
    return index_;
}


inline bool Foam::phasePair::const_iterator::operator==
(
    const const_iterator& iter
) const
{
    return (this->index_ == iter.index_);
}


inline bool Foam::phasePair::const_iterator::operator!=
(
    const const_iterator& iter
) const
{
    return !(this->operator==(iter));
}


inline const Foam::phaseModel&
Foam::phasePair::const_iterator::operator*() const
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
Foam::phasePair::const_iterator::operator()() const
{
    return operator*();
}


inline const Foam::phaseModel&
Foam::phasePair::const_iterator::otherPhase() const
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


inline Foam::phasePair::const_iterator&
Foam::phasePair::const_iterator::operator++()
{
    index_++;
    return *this;
}


inline Foam::phasePair::const_iterator
Foam::phasePair::const_iterator::operator++(int)
{
    const_iterator old = *this;
    this->operator++();
    return old;
}


inline Foam::phasePair::const_iterator Foam::phasePair::cbegin() const
{
    return const_iterator(*this);
}


inline Foam::phasePair::const_iterator Foam::phasePair::cend() const
{
    return const_iterator(*this, 2);
}


inline Foam::phasePair::const_iterator Foam::phasePair::begin() const
{
    return const_iterator(*this);
}


inline Foam::phasePair::const_iterator Foam::phasePair::end() const
{
    return const_iterator(*this, 2);
}


// ************************************************************************* //
