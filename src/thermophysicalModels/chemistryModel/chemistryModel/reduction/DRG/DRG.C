/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

template<class ThermoType>
Foam::chemistryReductionMethods::DRG<ThermoType>::DRG
(
    const dictionary& dict,
    chemistryModel<ThermoType>& chemistry
)
:
    chemistryReductionMethod<ThermoType>(dict, chemistry),
    searchInitSet_()
{
    const wordHashSet initSet(this->coeffDict(dict).lookup("initialSet"));
    forAllConstIter(wordHashSet, initSet, iter)
    {
        searchInitSet_.append(chemistry.thermo().species()[iter.key()]);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::chemistryReductionMethods::DRG<ThermoType>::~DRG()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::chemistryReductionMethods::DRG<ThermoType>::reduceMechanism
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    List<label>& ctos,
    DynamicList<label>& stoc,
    const label li
)
{
    chemistryReductionMethod<ThermoType>::initReduceMechanism();

    scalarField c1(this->nSpecie()+2, 0.0);

    for(label i=0; i<this->nSpecie(); i++)
    {
        c1[i] = c[i];
    }

    c1[this->nSpecie()] = T;
    c1[this->nSpecie()+1] = p;

    // Compute the rAB matrix
    RectangularMatrix<scalar> rABNum(this->nSpecie(),this->nSpecie(),0.0);
    scalarField rABDen(this->nSpecie(),0.0);

    // Number of initialised rAB for each lines
    Field<label> NbrABInit(this->nSpecie(),0);

    // Position of the initialised rAB, -1 when not initialised
    RectangularMatrix<label> rABPos(this->nSpecie(), this->nSpecie(), -1);

    // Index of the other species involved in the rABNum
    RectangularMatrix<label> rABOtherSpec(this->nSpecie(), this->nSpecie(), -1);

    forAll(this->chemistry_.reactions(), i)
    {
        const Reaction<ThermoType>& R = this->chemistry_.reactions()[i];

        // For each reaction compute omegai
        scalar omegaf, omegar;
        const scalar omegai = R.omega(p, T, c1, li, omegaf, omegar);

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

            List<bool> deltaBi(this->nSpecie(), false);
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

                    // Test if this rAB is not initialised
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

    // Set all species to inactive and activate them according
    // to rAB and initial set
    for (label i=0; i<this->nSpecie(); i++)
    {
        this->activeSpecies_[i] = false;
    }

    FIFOStack<label> Q;

    // Initialise the list of active species with the search initiating set
    // (SIS)
    for (label i=0; i<searchInitSet_.size(); i++)
    {
        label q = searchInitSet_[i];
        this->activeSpecies_[q] = true;
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
                }
            }
        }
    }

    chemistryReductionMethod<ThermoType>::endReduceMechanism(ctos, stoc);
}


// ************************************************************************* //
