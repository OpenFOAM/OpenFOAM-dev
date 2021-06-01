/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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

#include "DRGEP.H"
#include "SortableListDRGEP.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::chemistryReductionMethods::DRGEP<ThermoType>::DRGEP
(
    const IOdictionary& dict,
    TDACChemistryModel<ThermoType>& chemistry
)
:
    chemistryReductionMethod<ThermoType>(dict, chemistry),
    searchInitSet_(),
    sC_(this->nSpecie_,0),
    sH_(this->nSpecie_,0),
    sO_(this->nSpecie_,0),
    sN_(this->nSpecie_,0),
    NGroupBased_(50)
{
    const wordHashSet initSet(this->coeffsDict_.lookup("initialSet"));
    forAllConstIter(wordHashSet, initSet, iter)
    {
        searchInitSet_.append(chemistry.mixture().species()[iter.key()]);
    }

    if (this->coeffsDict_.found("NGroupBased"))
    {
        NGroupBased_ = this->coeffsDict_.template lookup<label>("NGroupBased");
    }

    for (label i=0; i<this->nSpecie_; i++)
    {
        const List<specieElement>& curSpecieComposition =
            chemistry.mixture().specieComposition(i);

        // for all elements in the current species
        forAll(curSpecieComposition, j)
        {
            const specieElement& curElement =
                curSpecieComposition[j];
            if (curElement.name() == "C")
            {
                sC_[i] = curElement.nAtoms();
            }
            else if (curElement.name() == "H")
            {
                sH_[i] = curElement.nAtoms();
            }
            else if (curElement.name() == "O")
            {
                sO_[i] = curElement.nAtoms();
            }
            else if (curElement.name() == "N")
            {
                sN_[i] = curElement.nAtoms();
            }
            else
            {
                Info<< "element not considered"<< endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::chemistryReductionMethods::DRGEP<ThermoType>::~DRGEP()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::chemistryReductionMethods::DRGEP<ThermoType>::reduceMechanism
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
)
{
    scalarField& completeC(this->chemistry_.completeC());
    scalarField c1(this->chemistry_.nEqns(), 0.0);

    for (label i=0; i<this->nSpecie_; i++)
    {
        c1[i] = c[i];
        completeC[i] = c[i];
    }

    c1[this->nSpecie_] = T;
    c1[this->nSpecie_+1] = p;

    // Compute the rAB matrix
    RectangularMatrix<scalar> rABNum(this->nSpecie_,this->nSpecie_,0.0);
    scalarField PA(this->nSpecie_,0.0);
    scalarField CA(this->nSpecie_,0.0);

    // Number of initialised rAB for each lines
    Field<label> NbrABInit(this->nSpecie_,0);
    // Position of the initialised rAB, -1 when not initialised
    RectangularMatrix<label> rABPos(this->nSpecie_, this->nSpecie_, -1);
    // Index of the other species involved in the rABNum
    RectangularMatrix<label> rABOtherSpec(this->nSpecie_, this->nSpecie_, -1);

    scalarField omegaV(this->chemistry_.reactions().size());
    forAll(this->chemistry_.reactions(), i)
    {
        const Reaction<ThermoType>& R = this->chemistry_.reactions()[i];

        // for each reaction compute omegai
        scalar omegaf, omegar;
        const scalar omegai = R.omega(p, T, c1, li, omegaf, omegar);
        omegaV[i] = omegai;

        // then for each pair of species composing this reaction,
        // compute the rAB matrix (separate the numerator and
        // denominator)

        DynamicList<scalar> wA(R.lhs().size()+R.rhs().size());
        DynamicList<label> wAID(R.lhs().size()+R.rhs().size());
        forAll(R.lhs(), s)// compute rAB for all species in the left hand side
        {
            label ss = R.lhs()[s].index;
            scalar sl = -R.lhs()[s].stoichCoeff; // vAi = v''-v' => here -v'
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
            deltaBi[ss] = false;

            while(!usedIndex.empty())
            {
                label curIndex = usedIndex.pop();
                if (deltaBi[curIndex])
                {
                    // disable to avoid counting it more than once
                    deltaBi[curIndex] = false;
                    // test if this rAB is not initialised
                    if (rABPos(ss, curIndex)==-1)
                    {
                        rABPos(ss, curIndex) = NbrABInit[ss];
                        NbrABInit[ss]++;
                        rABNum(ss, rABPos(ss, curIndex)) = sl*omegai;
                        rABOtherSpec(ss, rABPos(ss, curIndex)) = curIndex;
                    }
                    else
                    {
                        rABNum(ss, rABPos(ss, curIndex)) += sl*omegai;
                    }
                }
            }
            bool found(false);
            forAll(wAID, id)
            {
                if (ss==wAID[id])
                {
                    wA[id] += sl*omegai;
                    found = true;
                }
            }
            if (!found)
            {
                wA.append(sl*omegai);
                wAID.append(ss);
            }
        }

        // Compute rAB for all species in the right hand side
        forAll(R.rhs(), s)
        {
            label ss = R.rhs()[s].index;
            scalar sl = R.rhs()[s].stoichCoeff; // vAi = v''-v' => here v''
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
            deltaBi[ss] = false;

            while(!usedIndex.empty())
            {
                label curIndex = usedIndex.pop();
                if (deltaBi[curIndex])
                {
                    // disable to avoid counting it more than once
                    deltaBi[curIndex] = false;
                    // test if this rAB is not initialised
                    if (rABPos(ss, curIndex)==-1)
                    {
                        rABPos(ss, curIndex) = NbrABInit[ss];
                        NbrABInit[ss]++;
                        rABNum(ss, rABPos(ss, curIndex)) = sl*omegai;
                        rABOtherSpec(ss, rABPos(ss, curIndex)) = curIndex;
                    }
                    else
                    {
                        rABNum(ss, rABPos(ss, curIndex)) += sl*omegai;
                    }
                }
            }

            bool found(false);
            forAll(wAID, id)
            {
                if (ss==wAID[id])
                {
                    wA[id] += sl*omegai;
                    found = true;
                }
            }
            if (!found)
            {
                wA.append(sl*omegai);
                wAID.append(ss);
            }
        }

        wAID.shrink();
        // Now that every species of the reactions has been visited, we can
        // compute the production and consumption rate. This way, it avoids
        // getting wrong results when species are present in both lhs and rhs
        forAll(wAID, id)
        {
            if (wA[id] > 0.0)
            {
                if (PA[wAID[id]] == 0.0)
                {
                    PA[wAID[id]] = wA[id];
                }
                else
                {
                    PA[wAID[id]] += wA[id];
                }
            }
            else
            {
                if (CA[wAID[id]] == 0.0)
                {
                    CA[wAID[id]] = -wA[id];
                }
                else
                {
                    CA[wAID[id]] += -wA[id];
                }
            }
        }
    }
    // rii = 0.0 by definition

    // Compute the production rate of each element Pa
    label nElements = 4; // 4 main elements (C, H, O, N)
    scalarList Pa(nElements,0.0);
    scalarList Ca(nElements,0.0);

    // for (label q=0; q<SIS.size(); q++)
    for (label i=0; i<this->nSpecie_; i++)
    {
        Pa[0] += sC_[i]*max(0.0,(PA[i]-CA[i]));
        Ca[0] += sC_[i]*max(0.0,-(PA[i]-CA[i]));
        Pa[1] += sH_[i]*max(0.0,(PA[i]-CA[i]));
        Ca[1] += sH_[i]*max(0.0,-(PA[i]-CA[i]));
        Pa[2] += sO_[i]*max(0.0,(PA[i]-CA[i]));
        Ca[2] += sO_[i]*max(0.0,-(PA[i]-CA[i]));
        Pa[3] += sN_[i]*max(0.0,(PA[i]-CA[i]));
        Ca[3] += sN_[i]*max(0.0,-(PA[i]-CA[i]));
    }

    // Using the rAB matrix (numerator and denominator separated)
    // compute the R value according to the search initiating set
    scalarField Rvalue(this->nSpecie_,0.0);
    label speciesNumber = 0;
    List<bool> disabledSpecies(this->nSpecie_,false);

    // set all species to inactive and activate them according
    // to rAB and initial set
    for (label i=0; i<this->nSpecie_; i++)
    {
        this->activeSpecies_[i] = false;
    }
    // Initialise the FIFOStack for search set
    FIFOStack<label> Q;
    const labelList& SIS(this->searchInitSet_);
    DynamicList<label> QStart(SIS.size());
    DynamicList<scalar> alphaQ(SIS.size());

    // Compute the alpha coefficient and initialise the R value of the species
    // in the SIS
    for (label i=0; i<SIS.size(); i++)
    {
        label q = SIS[i];
        // compute alpha
        scalar alphaA(0.0);
        // C
        if (Pa[0] > vSmall)
        {
            scalar alphaTmp = (sC_[q]*mag(PA[q]-CA[q])/Pa[0]);
            if (alphaTmp > alphaA)
            {
                alphaA = alphaTmp;
            }
        }
        // H
        if (Pa[1] > vSmall)
        {
            scalar alphaTmp = (sH_[q]*mag(PA[q]-CA[q])/Pa[1]);
            if (alphaTmp > alphaA)
            {
                alphaA = alphaTmp;
            }
        }
        // O
        if (Pa[2] > vSmall)
        {
            scalar alphaTmp = (sO_[q]*mag(PA[q]-CA[q])/Pa[2]);
            if (alphaTmp > alphaA)
            {
                alphaA = alphaTmp;
            }
        }
        // N
        if (Pa[3] > vSmall)
        {
            scalar alphaTmp = (sN_[q]*mag(PA[q]-CA[q])/Pa[3]);
            if (alphaTmp > alphaA)
            {
                alphaA = alphaTmp;
            }
        }
        if (alphaA > this->tolerance())
        {
            this->activeSpecies_[q] = true;
            speciesNumber++;
            Q.push(q);
            QStart.append(q);
            alphaQ.append(1.0);
            Rvalue[q] = 1.0;
        }
        else
        {
            Rvalue[q] = alphaA;
        }
    }

    // if all species from the SIS has been removed
    // force the use of the species with maximum Rvalue
    if (Q.empty())
    {
        scalar Rmax=0.0;
        label specID=-1;
        forAll(SIS, i)
        {
            if (Rvalue[SIS[i]] > Rmax)
            {
                Rmax = Rvalue[SIS[i]];
                specID=SIS[i];
            }
        }
        Q.push(specID);
        QStart.append(specID);
        alphaQ.append(1.0);
        speciesNumber++;
        Rvalue[specID] = 1.0;
        this->activeSpecies_[specID] = true;
    }

    // Execute the main loop for R-value
    while (!Q.empty())
    {
        label u = Q.pop();
        scalar Den = max(PA[u],CA[u]);
        if (Den > vSmall)
        {
            for (label v=0; v<NbrABInit[u]; v++)
            {
                label otherSpec = rABOtherSpec(u, v);
                scalar rAB = mag(rABNum(u, v))/Den;
                if (rAB > 1)
                {
                    rAB = 1;
                }

                scalar Rtemp = Rvalue[u]*rAB;
                // a link analysed previously is stronger
                if (Rvalue[otherSpec] < Rtemp)
                {
                    Rvalue[otherSpec] = Rtemp;
                    // the (composed) link is stronger than the tolerance
                    if (Rtemp >= this->tolerance())
                    {
                        Q.push(otherSpec);
                        if (!this->activeSpecies_[otherSpec])
                        {
                            this->activeSpecies_[otherSpec] = true;
                            speciesNumber++;
                        }
                    }
                }
            }
        }
    }

    // Group-based reduction
    // number of species disabled in the first step
    label NDisabledSpecies(this->nSpecie_-speciesNumber);

    // while the number of removed species is greater than NGroupBased, the rAB
    // are reevaluated according to the group based definition for each loop the
    // temporary disabled species (in the first reduction) are sorted to disable
    // definitely the NGroupBased species with lower R then these R value a
    // reevaluated taking into account these disabled species
    while(NDisabledSpecies > NGroupBased_)
    {
        // keep track of disabled species using sortablelist to extract only
        // NGroupBased lower R value
        SortableListDRGEP<scalar> Rdisabled(NDisabledSpecies);
        labelList Rindex(NDisabledSpecies);
        label nD = 0;
        forAll(disabledSpecies, i)
        {
            // if just disabled and not in a previous loop
            if (!this->activeSpecies_[i] && !disabledSpecies[i])
            {
                // Note: non-reached species will be removed first (Rvalue=0)
                Rdisabled[nD] = Rvalue[i];
                Rindex[nD++] = i;
            }
        }
        // sort the Rvalue to obtain the NGroupBased lower R value
        Rdisabled.partialSort(NGroupBased_);
        labelList tmpIndex(Rdisabled.indices());

        // disable definitely NGroupBased species in this loop
        for (label i=0; i<NGroupBased_; i++)
        {
            disabledSpecies[Rindex[tmpIndex[i]]] = true;
        }
        NDisabledSpecies -= NGroupBased_;

        // reevaluate the rAB according to the group-based definition rAB{S} [1]
        // only update the numerator
        forAll(NbrABInit, i)
        {
            for (label v=0; v<NbrABInit[i]; v++)
            {
                rABNum(i, v) = 0.0;
            }
        }
        forAll(this->chemistry_.reactions(), i)
        {
            const Reaction<ThermoType>& R = this->chemistry_.reactions()[i];

            scalar omegai = omegaV[i];

            forAll(R.lhs(), s)
            {
                label ss = R.lhs()[s].index;
                scalar sl = -R.lhs()[s].stoichCoeff; // vAi = v''-v' => here -v'
                List<bool> deltaBi(this->nSpecie_, false);
                bool alreadyDisabled(false);
                FIFOStack<label> usedIndex;
                forAll(R.lhs(), j)
                {
                    label sj = R.lhs()[j].index;
                    usedIndex.push(sj);
                    deltaBi[sj] = true;
                    if (disabledSpecies[sj])
                    {
                        alreadyDisabled=true;
                    }
                }
                forAll(R.rhs(), j)
                {
                    label sj = R.rhs()[j].index;
                    usedIndex.push(sj);
                    deltaBi[sj] = true;
                    if (disabledSpecies[sj])
                    {
                        alreadyDisabled=true;
                    }
                }

                deltaBi[ss] = false;

                if (alreadyDisabled)
                {
                    // if one of the species in this reaction is disabled, all
                    // species connected to species ss are modified
                    for (label v=0; v<NbrABInit[ss]; v++)
                    {
                        rABNum(ss, v) += sl*omegai;
                    }
                }
                else
                {
                    while(!usedIndex.empty())
                    {
                        label curIndex = usedIndex.pop();
                        if (deltaBi[curIndex])
                        {
                            // disable to avoid counting it more than once
                            deltaBi[curIndex] = false;
                            rABNum(ss, rABPos(ss, curIndex)) += sl*omegai;
                        }
                    }
                }
            }

            forAll(R.rhs(), s)
            {
                label ss = R.rhs()[s].index;
                scalar sl = R.rhs()[s].stoichCoeff; // vAi = v''-v' => here v''
                List<bool> deltaBi(this->nSpecie_, false);
                bool alreadyDisabled(false);
                FIFOStack<label> usedIndex;
                forAll(R.lhs(), j)
                {
                    label sj = R.lhs()[j].index;
                    usedIndex.push(sj);
                    deltaBi[sj] = true;
                    if (disabledSpecies[sj])
                    {
                        alreadyDisabled=true;
                    }
                }
                forAll(R.rhs(), j)
                {
                    label sj = R.rhs()[j].index;
                    usedIndex.push(sj);
                    deltaBi[sj] = true;
                    if (disabledSpecies[sj])
                    {
                        alreadyDisabled=true;
                    }
                }

                deltaBi[ss] = false;

                if (alreadyDisabled)
                {
                    // if one of the species in this reaction is disabled, all
                    // species connected to species ss are modified
                    for (label v=0; v<NbrABInit[ss]; v++)
                    {
                        rABNum(ss, v) += sl*omegai;
                    }
                }
                else
                {
                    while(!usedIndex.empty())
                    {
                        label curIndex = usedIndex.pop();
                        if (deltaBi[curIndex])
                        {
                            deltaBi[curIndex] = false;
                            rABNum(ss, rABPos(ss, curIndex)) += sl*omegai;
                        }
                    }
                }
            }
        }

        forAll(QStart, qi)
        {
            label u = QStart[qi];
            Q.push(u);
        }

        while (!Q.empty())
        {
            label u = Q.pop();
            scalar Den = max(PA[u],CA[u]);
            if (Den!=0.0)
            {
                for (label v=0; v<NbrABInit[u]; v++)
                {
                    label otherSpec = rABOtherSpec(u, v);
                    if (!disabledSpecies[otherSpec])
                    {
                        scalar rAB = mag(rABNum(u, v))/Den;
                        if (rAB > 1)
                        {
                            rAB = 1;
                        }

                        scalar Rtemp = Rvalue[u]*rAB;
                        // a link analysed previously is stronger
                        if (Rvalue[otherSpec] < Rtemp)
                        {
                            Rvalue[otherSpec] = Rtemp;
                            if (Rtemp >= this->tolerance())
                            {
                                Q.push(otherSpec);
                                if (!this->activeSpecies_[otherSpec])
                                {
                                    this->activeSpecies_[otherSpec] = true;
                                    speciesNumber++;
                                    NDisabledSpecies--;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // End of group-based reduction

    // Put a flag on the reactions containing at least one removed species
    forAll(this->chemistry_.reactions(), i)
    {
        const Reaction<ThermoType>& R = this->chemistry_.reactions()[i];
        this->chemistry_.reactionsDisabled()[i] = false;
        forAll(R.lhs(), s)
        {
            label ss = R.lhs()[s].index;
            if (!this->activeSpecies_[ss])
            {
                this->chemistry_.reactionsDisabled()[i] = true;
                break;
            }
        }
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
    scalarField& simplifiedC(this->chemistry_.simplifiedC());
    simplifiedC.setSize(this->NsSimp_+2);
    DynamicList<label>& s2c(this->chemistry_.simplifiedToCompleteIndex());
    s2c.setSize(this->NsSimp_);
    Field<label>& c2s(this->chemistry_.completeToSimplifiedIndex());

    label j = 0;
    for (label i=0; i<this->nSpecie_; i++)
    {
        if (this->activeSpecies_[i])
        {
            s2c[j] = i;
            simplifiedC[j] = c[i];
            c2s[i] = j++;
            if (!this->chemistry_.active(i))
            {
                this->chemistry_.setActive(i);
            }
        }
        else
        {
            c2s[i] = -1;
        }
    }
    simplifiedC[this->NsSimp_] = T;
    simplifiedC[this->NsSimp_+1] = p;
    this->chemistry_.setNsDAC(this->NsSimp_);
    // change temporary Ns in chemistryModel
    // to make the function nEqns working
    this->chemistry_.setNSpecie(this->NsSimp_);
}


// ************************************************************************* //
