/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "chemistryReduction.H"
#include "chemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::chemistryReduction<ThermoType>::chemistryReduction
(
    const Foam::IOdictionary& dict,
    Foam::chemistryModel<ThermoType>& chemistry
)
:
    chemistryReductionMethod<ThermoType>(dict, chemistry),
    coeffsDict_(dict.subDict("reduction")),
    chemistry_(chemistry),
    activeSpecies_(chemistry.nSpecie(), false),
    log_(coeffsDict_.lookupOrDefault<Switch>("log", false)),
    tolerance_(coeffsDict_.lookupOrDefault<scalar>("tolerance", 1e-4)),
    clockTime_(clockTime()),
    sumnActiveSpecies_(0),
    sumn_(0),
    reduceMechCpuTime_(0)
{
    if (log_)
    {
        cpuReduceFile_ = chemistry.logFile("cpu_reduce.out");
        nActiveSpeciesFile_ = chemistry.logFile("nActiveSpecies.out");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::chemistryReduction<ThermoType>::~chemistryReduction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::chemistryReduction<ThermoType>::initReduceMechanism()
{
    if (log_)
    {
        clockTime_.timeIncrement();
    }
}


template<class ThermoType>
void Foam::chemistryReduction<ThermoType>::endReduceMechanism()
{
    // Change temporary Ns in chemistryModel
    // to make the function nEqns working
    this->chemistry_.setNSpecie(this->nActiveSpecies_);

    if (log_)
    {
        sumnActiveSpecies_ += this->nActiveSpecies_;
        sumn_++;
        reduceMechCpuTime_ += clockTime_.timeIncrement();
    }
}


template<class ThermoType>
void Foam::chemistryReduction<ThermoType>::update()
{
    if (log_)
    {
        cpuReduceFile_()
            << this->chemistry_.time().userTimeValue()
            << "    " << reduceMechCpuTime_ << endl;

        if (sumn_)
        {
            // Write average number of species
            nActiveSpeciesFile_()
                << this->chemistry_.time().userTimeValue()
                << "    " << sumnActiveSpecies_/sumn_ << endl;
        }

        sumnActiveSpecies_ = 0;
        sumn_ = 0;
        reduceMechCpuTime_ = 0;
    }
}


// ************************************************************************* //
