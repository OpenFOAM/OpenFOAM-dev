/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "multicomponentMixture.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::multicomponentMixture<ThermoType>::multicomponentMixture
(
    const dictionary& dict
)
:
    species_(dict.lookup("species")),
    specieThermos_(species_.size()),
    activeSpecies_(species_.size(), true)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::multicomponentMixture<ThermoType>::read(const dictionary& dict)
{
    forAll(species_, speciei)
    {
        specieThermos_.set
        (
            speciei,
            new ThermoType
            (
                species_[speciei],
                dict.subDict(species_[speciei])
            )
        );
    }
}


template<class ThermoType>
Foam::boolList& Foam::multicomponentMixture<ThermoType>::activeSpecies() const
{
    return activeSpecies_;
}


// ************************************************************************* //
