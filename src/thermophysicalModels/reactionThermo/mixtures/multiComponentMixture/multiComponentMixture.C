/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "multiComponentMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
Foam::PtrList<ThermoType>
Foam::multiComponentMixture<ThermoType>::readSpeciesData
(
    const dictionary& thermoDict
) const
{
    PtrList<ThermoType> specieThermos(species_.size());

    forAll(species_, i)
    {
        specieThermos.set
        (
            i,
            new ThermoType(thermoDict.subDict(species_[i]))
        );
    }

    return specieThermos;
}


template<class ThermoType>
typename Foam::multiComponentMixture<ThermoType>::speciesCompositionTable
Foam::multiComponentMixture<ThermoType>::readSpeciesComposition
(
    const dictionary& thermoDict,
    const speciesTable& species
) const
{
    speciesCompositionTable speciesComposition_;

    // Loop through all species in thermoDict to retrieve
    // the species composition
    forAll(species, si)
    {
        if (thermoDict.subDict(species[si]).isDict("elements"))
        {
            dictionary currentElements
            (
                thermoDict.subDict(species[si]).subDict("elements")
            );

            wordList currentElementsName(currentElements.toc());
            List<specieElement> currentComposition(currentElementsName.size());

            forAll(currentElementsName, eni)
            {
                currentComposition[eni].name() = currentElementsName[eni];

                currentComposition[eni].nAtoms() =
                    currentElements.lookupOrDefault
                    (
                        currentElementsName[eni],
                        0
                    );
            }

            // Add current specie composition to the hash table
            speciesCompositionTable::iterator specieCompositionIter
            (
                speciesComposition_.find(species[si])
            );

            if (specieCompositionIter != speciesComposition_.end())
            {
                speciesComposition_.erase(specieCompositionIter);
            }

            speciesComposition_.insert(species[si], currentComposition);
        }
    }

    return speciesComposition_;
}


template<class ThermoType>
void Foam::multiComponentMixture<ThermoType>::correctMassFractions()
{
    // Multiplication by 1.0 changes Yt patches to "calculated"
    volScalarField Yt("Yt", 1.0*Y_[0]);

    for (label n=1; n<Y_.size(); n++)
    {
        Yt += Y_[n];
    }

    if (mag(max(Yt).value()) < rootVSmall)
    {
        FatalErrorInFunction
            << "Sum of mass fractions is zero for species " << this->species()
            << exit(FatalError);
    }

    forAll(Y_, n)
    {
        Y_[n] /= Yt;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::multiComponentMixture<ThermoType>::multiComponentMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecieMixture
    (
        thermoDict,
        thermoDict.lookup("species"),
        mesh,
        phaseName
    ),
    specieThermos_(readSpeciesData(thermoDict)),
    speciesComposition_(readSpeciesComposition(thermoDict, species())),
    mixture_("mixture", specieThermos_[0]),
    mixtureVol_("volMixture", specieThermos_[0])
{
    correctMassFractions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::cellMixture
(
    const label celli
) const
{
    mixture_ = Y_[0][celli]*specieThermos_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixture_ += Y_[n][celli]*specieThermos_[n];
    }

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::patchFaceMixture
(
    const label patchi,
    const label facei
) const
{
    mixture_ = Y_[0].boundaryField()[patchi][facei]*specieThermos_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixture_ += Y_[n].boundaryField()[patchi][facei]*specieThermos_[n];
    }

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::cellVolMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    scalar rhoInv = 0.0;
    forAll(specieThermos_, i)
    {
        rhoInv += Y_[i][celli]/specieThermos_[i].rho(p, T);
    }

    mixtureVol_ =
        Y_[0][celli]/specieThermos_[0].rho(p, T)/rhoInv*specieThermos_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixtureVol_ +=
            Y_[n][celli]/specieThermos_[n].rho(p, T)/rhoInv*specieThermos_[n];
    }

    return mixtureVol_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::
patchFaceVolMixture
(
    const scalar p,
    const scalar T,
    const label patchi,
    const label facei
) const
{
    scalar rhoInv = 0.0;
    forAll(specieThermos_, i)
    {
        rhoInv +=
            Y_[i].boundaryField()[patchi][facei]/specieThermos_[i].rho(p, T);
    }

    mixtureVol_ =
        Y_[0].boundaryField()[patchi][facei]/specieThermos_[0].rho(p, T)/rhoInv
      * specieThermos_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixtureVol_ +=
            Y_[n].boundaryField()[patchi][facei]/specieThermos_[n].rho(p,T)
          / rhoInv*specieThermos_[n];
    }

    return mixtureVol_;
}


template<class ThermoType>
void Foam::multiComponentMixture<ThermoType>::read
(
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        specieThermos_[i] = ThermoType(thermoDict.subDict(species_[i]));
    }
}


// ************************************************************************* //
