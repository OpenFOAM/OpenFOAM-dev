/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "DRG.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryReductionMethods::DRG<CompType, ThermoType>::DRG
(
    const IOdictionary& dict,
    TDACChemistryModel<CompType, ThermoType>& chemistry
)
:
    chemistryReductionMethod<CompType, ThermoType>(dict, chemistry),
    searchInitSet_(this->coeffsDict_.subDict("initialSet").size())
{
    label j=0;
    dictionary initSet = this->coeffsDict_.subDict("initialSet");
    for (label i=0; i<chemistry.nSpecie(); i++)
    {
        if (initSet.found(chemistry.Y()[i].member()))
        {
            searchInitSet_[j++] = i;
        }
    }
    if (j<searchInitSet_.size())
    {
        FatalErrorInFunction
            << searchInitSet_.size()-j
            << " species in the initial set is not in the mechanism "
            << initSet
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryReductionMethods::DRG<CompType, ThermoType>::~DRG()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
void Foam::chemistryReductionMethods::DRG<CompType, ThermoType>::reduceMechanism
(
    const scalarField &c,
    const scalar T,
    const scalar p
)
{
    scalarField c1(this->nSpecie_+2, 0.0);

    for(label i=0; i<this->nSpecie_; i++)
    {
        c1[i] = c[i];
    }

    c1[this->nSpecie_] = T;
    c1[this->nSpecie_+1] = p;

    // Compute the rAB matrix
    RectangularMatrix<scalar> rABNum(this->nSpecie_,this->nSpecie_,0.0);
    scalarField rABDen(this->nSpecie_,0.0);

    // Number of initialized rAB for each lines
    Field<label> NbrABInit(this->nSpecie_,0);

    // Position of the initialized rAB, -1 when not initialized
    RectangularMatrix<label> rABPos(this->nSpecie_, this->nSpecie_, -1);

    // Index of the other species involved in the rABNum
    RectangularMatrix<label> rABOtherSpec(this->nSpecie_, this->nSpecie_, -1);

    scalar pf, cf, pr, cr;
    label lRef, rRef;
    forAll(this->chemistry_.reactions(), i)
    {
        const Reaction<ThermoType>& R = this->chemistry_.reactions()[i];
        // For each reaction compute omegai
        scalar omegai = this->chemistry_.omega
        (
         R, c1, T, p, pf, cf, lRef, pr, cr, rRef
         );


        // Then for each pair of species composing this reaction,
        // compute the rAB matrix (separate the numerator and
        // denominator)
        DynamicList<scalar> wA(R.lhs().size()+R.rhs().size());
        DynamicList<label> wAID(R.lhs().size()+R.rhs().size());

        forAll(R.lhs(), s)
        {
            label ss = R.lhs()[s].index;
            scalar sl = -R.lhs()[s].stoichCoeff;
            bool found(false);
            forAll(wAID, id)
            {
                if (ss==wAID[id])
                {
                    wA[id] += sl*omegai;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                wA.append(sl*omegai);
                wAID.append(ss);
            }
        }
        forAll(R.rhs(), s)
        {
            label ss = R.rhs()[s].index;
            scalar sl = R.rhs()[s].stoichCoeff;
            bool found(false);
            forAll(wAID, id)
            {
                if (ss==wAID[id])
                {
                    wA[id] += sl*omegai;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                wA.append(sl*omegai);
                wAID.append(ss);
            }
        }

        // Now that all nuAi*wi are computed, without counting twice species
        // present in both rhs and lhs, we can update the denominator and
        // numerator for the rAB
        wAID.shrink();
        forAll(wAID, id)
        {
            label curID = wAID[id];

            // Absolute value of aggregated value
            scalar curwA = ((wA[id]>=0) ? wA[id] : -wA[id]);

            List<bool> deltaBi(this->nSpecie_, false);
            FIFOStack<label> usedIndex;
            forAll(R.lhs(), j)
            {
                label sj = R.lhs()[j].index;
                usedIndex.push(sj);
                deltaBi[sj] = true;
            }
            forAll(R.rhs(), j)
            {
                label sj = R.rhs()[j].index;
                usedIndex.push(sj);
                deltaBi[sj] = true;
            }

            // Disable for self reference (by definition rAA=0)
            deltaBi[curID] = false;
            while(!usedIndex.empty())
            {
                label curIndex = usedIndex.pop();

                if (deltaBi[curIndex])
                {
                    // Disable to avoid counting it more than once
                    deltaBi[curIndex] = false;

                    // Test if this rAB is not initialized
                    if (rABPos(curID, curIndex)==-1)
                    {
                        rABPos(curID, curIndex) = NbrABInit[curID];
                        NbrABInit[curID]++;
                        rABNum(curID, rABPos(curID, curIndex)) = curwA;
                        rABOtherSpec(curID, rABPos(curID, curIndex)) = curIndex;
                    }
                    else
                    {
                        rABNum(curID, rABPos(curID, curIndex)) += curwA;
                    }
                }
            }
            if (rABDen[curID] == 0.0)
            {
                rABDen[curID] = curwA;
            }
            else
            {
                rABDen[curID] +=curwA;
            }
        }
    }
    // rii = 0.0 by definition

    label speciesNumber = 0;

    // Set all species to inactive and activate them according
    // to rAB and initial set
    for (label i=0; i<this->nSpecie_; i++)
    {
        this->activeSpecies_[i] = false;
    }

    FIFOStack<label> Q;

    // Initialize the list of active species with the search initiating set
    // (SIS)
    for (label i=0; i<searchInitSet_.size(); i++)
    {
        label q = searchInitSet_[i];
        this->activeSpecies_[q] = true;
        speciesNumber++;
        Q.push(q);
    }

    // Breadth first search with rAB
    while (!Q.empty())
    {
        label u = Q.pop();
        scalar Den = rABDen[u];

        if (Den > vSmall)
        {
            for (label v=0; v<NbrABInit[u]; v++)
            {
                label otherSpec = rABOtherSpec(u, v);
                scalar rAB = rABNum(u, v)/Den;

                if (rAB > 1)
                {
                    rAB = 1;
                }

                // Include B only if rAB is above the tolerance and if the
                // species was not searched before
                if
                (
                    rAB >= this->tolerance()
                 && !this->activeSpecies_[otherSpec]
                )
                {
                    Q.push(otherSpec);
                    this->activeSpecies_[otherSpec] = true;
                    speciesNumber++;
                }
            }
        }
    }

    // Put a flag on the reactions containing at least one removed species
    forAll(this->chemistry_.reactions(), i)
    {
        const Reaction<ThermoType>& R = this->chemistry_.reactions()[i];
        this->chemistry_.reactionsDisabled()[i] = false;

        forAll(R.lhs(), s)
        {
            label ss = R.lhs()[s].index;

            // The species is inactive then the reaction is removed
            if (!this->activeSpecies_[ss])
            {
                // Flag the reaction to disable it
                this->chemistry_.reactionsDisabled()[i] = true;
                break;
            }
        }

        // If the reaction has not been disabled yet
        if (!this->chemistry_.reactionsDisabled()[i])
        {
            forAll(R.rhs(), s)
            {
                label ss = R.rhs()[s].index;
                if (!this->activeSpecies_[ss])
                {
                    this->chemistry_.reactionsDisabled()[i] = true;
                    break;
                }
            }
        }
    }

    this->NsSimp_ = speciesNumber;
    this->chemistry_.simplifiedC().setSize(this->NsSimp_+2);
    this->chemistry_.simplifiedToCompleteIndex().setSize(this->NsSimp_);

    label j = 0;
    for (label i=0; i<this->nSpecie_; i++)
    {
        if (this->activeSpecies_[i])
        {
            this->chemistry_.simplifiedToCompleteIndex()[j] = i;
            this->chemistry_.simplifiedC()[j] = c[i];
            this->chemistry_.completeToSimplifiedIndex()[i] = j++;
            if (!this->chemistry_.active(i))
            {
                this->chemistry_.setActive(i);
            }
        }
        else
        {
            this->chemistry_.completeToSimplifiedIndex()[i] = -1;
        }
    }

    this->chemistry_.simplifiedC()[this->NsSimp_] = T;
    this->chemistry_.simplifiedC()[this->NsSimp_+1] = p;
    this->chemistry_.setNsDAC(this->NsSimp_);

    // Change temporary Ns in chemistryModel
    // to make the function nEqns working
    this->chemistry_.setNSpecie(this->NsSimp_);
}


// ************************************************************************* //
