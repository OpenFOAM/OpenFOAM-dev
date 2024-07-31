/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "chemistryReductionMethod.H"
#include "chemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::chemistryReductionMethod<ThermoType>::chemistryReductionMethod
(
    Foam::chemistryModel<ThermoType>& chemistry
)
:
    coeffDict_(),
    chemistry_(chemistry),
    nSpecie_(chemistry.nSpecie()),
    nActiveSpecies_(chemistry.nSpecie()),
    reactionsDisabled_(chemistry.nReaction(), false),
    activeSpecies_(chemistry.nSpecie(), true),
    log_(false),
    tolerance_(NaN),
    sumnActiveSpecies_(0),
    sumn_(0),
    reduceMechCpuTime_(0)
{}


template<class ThermoType>
Foam::chemistryReductionMethod<ThermoType>::chemistryReductionMethod
(
    const Foam::IOdictionary& dict,
    Foam::chemistryModel<ThermoType>& chemistry
)
:
    coeffDict_(dict.subDict("reduction")),
    chemistry_(chemistry),
    nSpecie_(chemistry.nSpecie()),
    nActiveSpecies_(chemistry.nSpecie()),
    reactionsDisabled_(chemistry.nReaction(), false),
    activeSpecies_(chemistry.nSpecie(), false),
    log_(coeffDict_.lookupOrDefault<Switch>("log", false)),
    tolerance_(coeffDict_.lookupOrDefault<scalar>("tolerance", 1e-4)),
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
Foam::chemistryReductionMethod<ThermoType>::~chemistryReductionMethod()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::chemistryReductionMethod<ThermoType>::initReduceMechanism()
{
    if (log_)
    {
        cpuTime_.cpuTimeIncrement();
    }
}


template<class ThermoType>
void Foam::chemistryReductionMethod<ThermoType>::endReduceMechanism
(
    List<label>& ctos,
    DynamicList<label>& stoc
)
{
    // Disable reactions containing removed species
    forAll(chemistry_.reactions(), i)
    {
        const Reaction<ThermoType>& R = chemistry_.reactions()[i];
        reactionsDisabled_[i] = false;

        forAll(R.lhs(), s)
        {
            label ss = R.lhs()[s].index;
            if (!activeSpecies_[ss])
            {
                reactionsDisabled_[i] = true;
                break;
            }
        }

        if (!reactionsDisabled_[i])
        {
            forAll(R.rhs(), s)
            {
                label ss = R.rhs()[s].index;
                if (!activeSpecies_[ss])
                {
                    reactionsDisabled_[i] = true;
                    break;
                }
            }
        }
    }

    // Set the total number of active species
    nActiveSpecies_ = count(activeSpecies_, true);

    // Set the indexing arrays
    stoc.setSize(nActiveSpecies_);
    for (label i=0, j=0; i<nSpecie(); i++)
    {
        if (activeSpecies_[i])
        {
            stoc[j] = i;
            ctos[i] = j++;
            if (!chemistry_.thermo().speciesActive()[i])
            {
                chemistry_.thermo().setSpecieActive(i);
            }
        }
        else
        {
            ctos[i] = -1;
        }
    }

    // Change the number of species in the chemistry model
    chemistry_.setNSpecie(nActiveSpecies_);

    if (log_)
    {
        sumnActiveSpecies_ += nActiveSpecies_;
        sumn_++;
        reduceMechCpuTime_ += cpuTime_.cpuTimeIncrement();
    }
}


template<class ThermoType>
void Foam::chemistryReductionMethod<ThermoType>::update()
{
    if (log_)
    {
        cpuReduceFile_()
            << chemistry_.time().userTimeValue()
            << "    " << reduceMechCpuTime_ << endl;

        if (sumn_)
        {
            // Write average number of species
            nActiveSpeciesFile_()
                << chemistry_.time().userTimeValue()
                << "    " << sumnActiveSpecies_/sumn_ << endl;
        }

        sumnActiveSpecies_ = 0;
        sumn_ = 0;
        reduceMechCpuTime_ = 0;
    }
}


// ************************************************************************* //
