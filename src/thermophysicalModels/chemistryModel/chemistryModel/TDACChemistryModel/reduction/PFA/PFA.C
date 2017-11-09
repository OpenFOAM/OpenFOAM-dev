/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
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

#include "PFA.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryReductionMethods::PFA<CompType, ThermoType>::PFA
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
            << " species in the intial set is not in the mechanism "
            << initSet
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryReductionMethods::PFA<CompType, ThermoType>::~PFA()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
void Foam::chemistryReductionMethods::PFA<CompType, ThermoType>::reduceMechanism
(
    const scalarField &c,
    const scalar T,
    const scalar p
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
    RectangularMatrix<scalar> PAB(this->nSpecie_,this->nSpecie_,0.0);
    RectangularMatrix<scalar> CAB(this->nSpecie_,this->nSpecie_,0.0);
    scalarField PA(this->nSpecie_,0.0);
    scalarField CA(this->nSpecie_,0.0);

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
        // for each reaction compute omegai
        scalar omegai = this->chemistry_.omega
        (
            R, c1, T, p, pf, cf, lRef, pr, cr, rRef
        );

        // then for each pair of species composing this reaction,
        // compute the rAB matrix (separate the numerator and
        // denominator)

        DynamicList<scalar> wA(R.lhs().size()+R.rhs().size());
        DynamicList<label> wAID(R.lhs().size()+R.rhs().size());

        forAll(R.lhs(), s)// compute rAB for all species in the left hand side
        {
            label ss = R.lhs()[s].index;
            scalar sl = -R.lhs()[s].stoichCoeff; // vAi = v''-v' => here -v'
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
        forAll(R.rhs(), s) // compute rAB for all species in the right hand side
        {
            label ss = R.rhs()[s].index;
            scalar sl = R.rhs()[s].stoichCoeff; // vAi = v''-v' => here v''
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

        wAID.shrink();
        forAll(wAID, id)
        {
            label curID = wAID[id];
            scalar curwA = wA[id];
            List<bool> deltaBi(this->nSpecie_, false);
            FIFOStack<label> usedIndex;
            forAll(R.lhs(),j)
            {
                label sj = R.lhs()[j].index;
                usedIndex.push(sj);
                deltaBi[sj] = true;
            }
            forAll(R.rhs(),j)
            {
                label sj = R.rhs()[j].index;
                usedIndex.push(sj);
                deltaBi[sj] = true;
            }

            deltaBi[curID] = false;

            while(!usedIndex.empty())
            {
                label curIndex = usedIndex.pop();

                if (deltaBi[curIndex])
                {
                    deltaBi[curIndex] = false;
                    if (rABPos(curID, curIndex)==-1)
                    {
                        rABPos(curID, curIndex) = NbrABInit[curID];
                        rABOtherSpec(curID, NbrABInit[curID]) = curIndex;
                        if (curwA > 0.0)
                        {
                            PAB(curID, NbrABInit[curID]) = curwA;
                        }
                        else
                        {
                            CAB(curID, NbrABInit[curID]) = -curwA;
                        }
                        NbrABInit[curID]++;
                    }
                    else
                    {
                        if (curwA > 0.0)
                        {
                            PAB(curID, rABPos(curID, curIndex)) += curwA;
                        }
                        else
                        {
                            CAB(curID, rABPos(curID, curIndex)) += -curwA;
                        }
                    }
                }
            }
            // Now that every species of the reactions has been visited, we can
            // compute the production and consumption rate. It avoids getting
            // wrong results when species are present in both lhs and rhs
            if (curwA > 0.0)
            {
                if (PA[curID] == 0.0)
                {
                    PA[curID] = curwA;
                }
                else
                {
                    PA[curID] += curwA;
                }
            }
            else
            {
                if (CA[curID] == 0.0)
                {
                    CA[curID] = -curwA;
                }
                else
                {
                    CA[curID] += -curwA;
                }
            }
        }
    }
    // rii = 0.0 by definition

    // Compute second generation link strength
    // For all species A look at all rAri of the connected species ri and
    // compute rriB with all the connected species of ri, B different of A.  If
    // a new species is connected, add it to the list of connected species.  It
    // is a connection of second generation and it will be aggregated in the
    // final step to evaluate the total connection strength (or path flux).
    // Compute rsecond=rAri*rriB with A!=ri!=B
    RectangularMatrix<scalar> PAB2nd(this->nSpecie_,this->nSpecie_,0.0);
    RectangularMatrix<scalar> CAB2nd(this->nSpecie_,this->nSpecie_,0.0);

    // Number of initialized rAB for each lines
    Field<label> NbrABInit2nd(this->nSpecie_, 0);

    // Position of the initialized rAB, -1 when not initialized
    RectangularMatrix<label> rABPos2nd(this->nSpecie_, this->nSpecie_, -1);

    // Index of the other species involved in the rABNum
    RectangularMatrix<label> rABOtherSpec2nd
    (
        this->nSpecie_, this->nSpecie_, -1
    );

    forAll(NbrABInit, A)
    {
        for (int i=0; i<NbrABInit[A]; i++)
        {
            label ri = rABOtherSpec(A, i);
            scalar maxPACA = max(PA[ri],CA[ri]);
            if (maxPACA > VSMALL)
            {
                for (int j=0; j<NbrABInit[ri]; j++)
                {
                    label B = rABOtherSpec(ri, j);
                    if (B != A) // if B!=A and also !=ri by definition
                    {
                        if (rABPos2nd(A, B)==-1)
                        {
                            rABPos2nd(A, B) = NbrABInit2nd[A]++;
                            rABOtherSpec2nd(A, rABPos2nd(A, B)) = B;
                            PAB2nd(A, rABPos2nd(A, B)) =
                                PAB(A, i)*PAB(ri, j)/maxPACA;
                            CAB2nd(A, rABPos2nd(A, B)) =
                                CAB(A, i)*CAB(ri, j)/maxPACA;
                        }
                        else
                        {
                            PAB2nd(A, rABPos2nd(A, B)) +=
                                PAB(A, i)*PAB(ri, j)/maxPACA;
                            CAB2nd(A, rABPos2nd(A, B)) +=
                                CAB(A, i)*CAB(ri, j)/maxPACA;
                        }
                    }
                }
            }
        }
    }

    // Using the rAB matrix (numerator and denominator separated)
    label speciesNumber = 0;

    // set all species to inactive and activate them according
    // to rAB and initial set
    for (label i=0; i<this->nSpecie_; i++)
    {
        this->activeSpecies_[i] = false;
    }

    // Initialize the FIFOStack for search set
    const labelList& SIS(this->searchInitSet_);
    FIFOStack<label> Q;

    for (label i=0; i<SIS.size(); i++)
    {
        label q = SIS[i];
        this->activeSpecies_[q] = true;
        speciesNumber++;
        Q.push(q);
    }

    // Execute the main loop for R-value
    while (!Q.empty())
    {
        label u = Q.pop();
        scalar Den = max(PA[u],CA[u]);

        if (Den!=0.0)
        {
            // first generation
            for (label v=0; v<NbrABInit[u]; v++)
            {
                label otherSpec = rABOtherSpec(u, v);
                scalar rAB = (PAB(u, v)+CAB(u, v))/Den; // first generation
                label id2nd = rABPos2nd(u, otherSpec);
                if (id2nd !=-1)// if there is a second generation link
                {
                    rAB += (PAB2nd(u, id2nd)+CAB2nd(u, id2nd))/Den;
                }
                // the link is stronger than the user-defined tolerance
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
            // second generation link only (for those without first link)
            for (label v=0; v<NbrABInit2nd[u]; v++)
            {
                label otherSpec = rABOtherSpec2nd(u, v);
                scalar rAB = (PAB2nd(u, v)+CAB2nd(u, v))/Den;
                // the link is stronger than the user-defined tolerance
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
