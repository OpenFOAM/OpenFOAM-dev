/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2026 OpenFOAM Foundation
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

#include "PhaseSurface_reactionRate.H"
#include "diameterModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Rate>
inline Foam::reactionRates::PhaseSurface<Rate>::PhaseSurface
(
    const speciesTable& species,
    const objectRegistry& ob,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    k_(species, dims*dimensions::length, dict),
    phaseName_(dict.lookup("phase")),
    ob_(ob),
    tAv_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Rate>
inline void Foam::reactionRates::PhaseSurface<Rate>::preEvaluate() const
{
    k_.preEvaluate();

    const phaseModel& phase =
        ob_.lookupObject<phaseModel>(IOobject::groupName("alpha", phaseName_));

    tAv_ = phase.diameter().Av();
}


template<class Rate>
inline void Foam::reactionRates::PhaseSurface<Rate>::postEvaluate() const
{
    k_.postEvaluate();

    tAv_.clear();
}


template<class Rate>
inline Foam::scalar Foam::reactionRates::PhaseSurface<Rate>::k
(
    const reaction::state& s
) const
{
    return k_.k(s)*tAv_()[s.li];
}


template<class Rate>
inline Foam::scalar Foam::reactionRates::PhaseSurface<Rate>::ddT
(
    const reaction::state& s
) const
{
    return k_.ddT(s)*tAv_()[s.li];
}


template<class Rate>
inline void Foam::reactionRates::PhaseSurface<Rate>::ddc
(
    const reaction::state& s,
    scalarField& ddc
) const
{
    k_.ddc(s, ddc);
    ddc *= tAv_()[s.li];
}


template<class Rate>
inline void Foam::reactionRates::PhaseSurface<Rate>::write(Ostream& os) const
{
    k_.write(os);
    writeEntry(os, "phase", phaseName_);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Rate>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const reactionRates::PhaseSurface<Rate>& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
