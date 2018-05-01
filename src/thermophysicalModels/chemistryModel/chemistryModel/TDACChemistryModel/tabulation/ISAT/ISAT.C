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

#include "ISAT.H"
#include "LUscalarMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::ISAT
(
    const dictionary& chemistryProperties,
    TDACChemistryModel<CompType, ThermoType>& chemistry
)
:
    chemistryTabulationMethod<CompType, ThermoType>
    (
        chemistryProperties,
        chemistry
    ),
    chemisTree_(chemistry, this->coeffsDict_),
    scaleFactor_(chemistry.nEqns() + ((this->variableTimeStep()) ? 1 : 0), 1),
    runTime_(chemistry.time()),
    chPMaxLifeTime_
    (
        this->coeffsDict_.lookupOrDefault("chPMaxLifeTime", INT_MAX)
    ),
    maxGrowth_(this->coeffsDict_.lookupOrDefault("maxGrowth", INT_MAX)),
    checkEntireTreeInterval_
    (
        this->coeffsDict_.lookupOrDefault("checkEntireTreeInterval", INT_MAX)
    ),
    maxDepthFactor_
    (
        this->coeffsDict_.lookupOrDefault
        (
            "maxDepthFactor",
            (chemisTree_.maxNLeafs() - 1)
           /(log(scalar(chemisTree_.maxNLeafs()))/log(2.0))
        )
    ),
    minBalanceThreshold_
    (
        this->coeffsDict_.lookupOrDefault
        (
            "minBalanceThreshold",0.1*chemisTree_.maxNLeafs()
        )
    ),
    MRURetrieve_(this->coeffsDict_.lookupOrDefault("MRURetrieve", false)),
    maxMRUSize_(this->coeffsDict_.lookupOrDefault("maxMRUSize", 0)),
    lastSearch_(nullptr),
    growPoints_(this->coeffsDict_.lookupOrDefault("growPoints", true)),
    nRetrieved_(0),
    nGrowth_(0),
    nAdd_(0),
    cleaningRequired_(false)
{
    if (this->active_)
    {
        dictionary scaleDict(this->coeffsDict_.subDict("scaleFactor"));
        label Ysize = this->chemistry_.Y().size();
        scalar otherScaleFactor = readScalar(scaleDict.lookup("otherSpecies"));
        for (label i=0; i<Ysize; i++)
        {
            if (!scaleDict.found(this->chemistry_.Y()[i].member()))
            {
                scaleFactor_[i] = otherScaleFactor;
            }
            else
            {
                scaleFactor_[i] =
                    readScalar
                    (
                        scaleDict.lookup(this->chemistry_.Y()[i].member())
                    );
            }
        }
        scaleFactor_[Ysize] = readScalar(scaleDict.lookup("Temperature"));
        scaleFactor_[Ysize + 1] = readScalar(scaleDict.lookup("Pressure"));
        if (this->variableTimeStep())
        {
            scaleFactor_[Ysize + 2] = readScalar(scaleDict.lookup("deltaT"));
        }
    }

    if (this->variableTimeStep())
    {
        nAdditionalEqns_ = 3;
    }
    else
    {
        nAdditionalEqns_ = 2;
    }

    if (this->log())
    {
        nRetrievedFile_ = chemistry.logFile("found_isat.out");
        nGrowthFile_ = chemistry.logFile("growth_isat.out");
        nAddFile_ = chemistry.logFile("add_isat.out");
        sizeFile_ = chemistry.logFile("size_isat.out");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::~ISAT()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CompType, class ThermoType>
void Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::addToMRU
(
    chemPointISAT<CompType, ThermoType>* phi0
)
{
    if (maxMRUSize_ > 0 && MRURetrieve_)
    {
        // First search if the chemPoint is already in the list
        bool isInList = false;
        typename SLList <chemPointISAT<CompType, ThermoType>*>::iterator iter =
            MRUList_.begin();
        for ( ; iter != MRUList_.end(); ++iter)
        {
            if (iter() == phi0)
            {
                isInList = true;
                break;
            }
        }
        // If it is in the list, then move it to front
        if (isInList)
        {
            if (iter() != MRUList_.first())
            {
                // iter hold the position of the element to move
                MRUList_.remove(iter);

                // Insert the element in front of the list
                MRUList_.insert(phi0);
            }
        }
        else // chemPoint not yet in the list, iter is last
        {
            if (MRUList_.size() == maxMRUSize_)
            {
                if (iter() == MRUList_.last())
                {
                    MRUList_.remove(iter);
                    MRUList_.insert(phi0);
                }
                else
                {
                    FatalErrorInFunction
                        << "wrong MRUList construction"
                        << exit(FatalError);
                }
            }
            else
            {
                MRUList_.insert(phi0);
            }
        }
    }
}


template<class CompType, class ThermoType>
void Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::calcNewC
(
    chemPointISAT<CompType, ThermoType>* phi0,
    const scalarField& phiq,
    scalarField& Rphiq
)
{
    label nEqns = this->chemistry_.nEqns(); // Species, T, p
    bool mechRedActive = this->chemistry_.mechRed()->active();
    Rphiq = phi0->Rphi();
    scalarField dphi(phiq-phi0->phi());
    const scalarSquareMatrix& gradientsMatrix = phi0->A();
    List<label>& completeToSimplified(phi0->completeToSimplifiedIndex());

    // Rphiq[i]=Rphi0[i]+A(i, j)dphi[j]
    // where Aij is dRi/dphi_j
    for (label i=0; i<nEqns-nAdditionalEqns_; i++)
    {
        if (mechRedActive)
        {
            label si = completeToSimplified[i];
            // The species is active
            if (si != -1)
            {
                for (label j=0; j<nEqns-2; j++)
                {
                    label sj = completeToSimplified[j];
                    if (sj != -1)
                    {
                        Rphiq[i] += gradientsMatrix(si, sj)*dphi[j];
                    }
                }
                Rphiq[i] +=
                    gradientsMatrix(si, phi0->nActiveSpecies())*dphi[nEqns - 2];
                Rphiq[i] +=
                    gradientsMatrix(si, phi0->nActiveSpecies() + 1)
                   *dphi[nEqns - 1];

                if (this->variableTimeStep())
                {
                    Rphiq[i] +=
                        gradientsMatrix(si, phi0->nActiveSpecies() + 2)
                       *dphi[nEqns];
                }

                // As we use an approximation of A, Rphiq should be checked for
                // negative values
                Rphiq[i] = max(0.0,Rphiq[i]);
            }
            // The species is not active A(i, j) = I(i, j)
            else
            {
                Rphiq[i] += dphi[i];
                Rphiq[i] = max(0.0,Rphiq[i]);
            }
        }
        else // Mechanism reduction is not active
        {
            for (label j=0; j<nEqns; j++)
            {
                Rphiq[i] += gradientsMatrix(i, j)*dphi[j];
            }
            // As we use a first order gradient matrix, Rphiq should be checked
            // for negative values
            Rphiq[i] = max(0.0,Rphiq[i]);
        }
    }
}


template<class CompType, class ThermoType>
bool Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::grow
(
    chemPointISAT<CompType, ThermoType>* phi0,
    const scalarField& phiq,
    const scalarField& Rphiq
)
{
    // If the pointer to the chemPoint is nullptr, the function stops
    if (!phi0)
    {
        return false;
    }

    // Raise a flag when the chemPoint used has been grown more than the
    // allowed number of time
    if (phi0->nGrowth() > maxGrowth_)
    {
        cleaningRequired_ = true;
        phi0->toRemove() = true;
        return false;
    }

    // If the solution RphiQ is still within the tolerance we try to grow it
    // in some cases this might result in a failure and the grow function of
    // the chemPoint returns false
    if (phi0->checkSolution(phiq, Rphiq))
    {
        return phi0->grow(phiq);
    }
    // The actual solution and the approximation given by ISAT are too different
    else
    {
        return false;
    }
}


template<class CompType, class ThermoType>
bool
Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::cleanAndBalance()
{
    bool treeModified(false);

    // Check all chemPoints to see if we need to delete some of the chemPoints
    // according to the ellapsed time and number of growths
    chemPointISAT<CompType, ThermoType>* x = chemisTree_.treeMin();
    while(x != nullptr)
    {
        chemPointISAT<CompType, ThermoType>* xtmp =
            chemisTree_.treeSuccessor(x);

        scalar elapsedTimeSteps = this->chemistry_.timeSteps() - x->timeTag();

        if ((elapsedTimeSteps > chPMaxLifeTime_) || (x->nGrowth() > maxGrowth_))
        {
            chemisTree_.deleteLeaf(x);
            treeModified = true;
        }
        x = xtmp;
    }
    // Check if the tree should be balanced according to criterion:
    //  -the depth of the tree bigger than a*log2(size), log2(size) being the
    //      ideal depth (e.g. 4 leafs can be stored in a tree of depth 2)
    if
    (
        chemisTree_.size() > minBalanceThreshold_
     && chemisTree_.depth() >
        maxDepthFactor_*log(scalar(chemisTree_.size()))/log(2.0)
    )
    {
        chemisTree_.balance();
        MRUList_.clear();
        treeModified = true;
    }

    // Return a bool to specify if the tree structure has been modified and is
    // now below the user specified limit (true if not full)
    return (treeModified && !chemisTree_.isFull());
}


template<class CompType, class ThermoType>
void Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::computeA
(
    scalarSquareMatrix& A,
    const scalarField& Rphiq,
    const scalar rhoi,
    const scalar dt
)
{
    bool mechRedActive = this->chemistry_.mechRed()->active();
    label speciesNumber = this->chemistry_.nSpecie();
    scalarField Rcq(this->chemistry_.nEqns() + nAdditionalEqns_ - 2);
    for (label i=0; i<speciesNumber; i++)
    {
        label s2c = i;
        if (mechRedActive)
        {
            s2c = this->chemistry_.simplifiedToCompleteIndex()[i];
        }
        Rcq[i] = rhoi*Rphiq[s2c]/this->chemistry_.specieThermo()[s2c].W();
    }
    Rcq[speciesNumber] = Rphiq[Rphiq.size() - nAdditionalEqns_];
    Rcq[speciesNumber+1] = Rphiq[Rphiq.size() - nAdditionalEqns_ + 1];
    if (this->variableTimeStep())
    {
        Rcq[speciesNumber + 2] = Rphiq[Rphiq.size() - nAdditionalEqns_ + 2];
    }

    // Aaa is computed implicitly,
    // A is given by A = C(psi0, t0+dt), where C is obtained through solving
    // d/dt C(psi0,t) = J(psi(t))C(psi0,t)
    // If we solve it implicitly:
    // (C(psi0, t0+dt) - C(psi0,t0))/dt = J(psi(t0+dt))C(psi0,t0+dt)
    // The Jacobian is thus computed according to the mapping
    // C(psi0,t0+dt)*(I-dt*J(psi(t0+dt))) = C(psi0, t0)
    // A = C(psi0,t0)/(I-dt*J(psi(t0+dt)))
    // where C(psi0,t0) = I

    this->chemistry_.jacobian(runTime_.value(), Rcq, A);

    // The jacobian is computed according to the molar concentration
    // the following conversion allows the code to use A with mass fraction
    for (label i=0; i<speciesNumber; i++)
    {
        label si = i;

        if (mechRedActive)
        {
            si = this->chemistry_.simplifiedToCompleteIndex()[i];
        }

        for (label j=0; j<speciesNumber; j++)
        {
            label sj = j;
            if (mechRedActive)
            {
                sj = this->chemistry_.simplifiedToCompleteIndex()[j];
            }
            A(i, j) *=
              -dt*this->chemistry_.specieThermo()[si].W()
               /this->chemistry_.specieThermo()[sj].W();
        }

        A(i, i) += 1;
        // Columns for pressure and temperature
        A(i, speciesNumber) *=
            -dt*this->chemistry_.specieThermo()[si].W()/rhoi;
        A(i, speciesNumber+1) *=
            -dt*this->chemistry_.specieThermo()[si].W()/rhoi;
    }

    // For temperature and pressure, only unity on the diagonal
    A(speciesNumber, speciesNumber) = 1;
    A(speciesNumber + 1, speciesNumber + 1) = 1;
    if (this->variableTimeStep())
    {
        A[speciesNumber + 2][speciesNumber + 2] = 1;
    }

    // Inverse of (I-dt*J(psi(t0+dt)))
    LUscalarMatrix LUA(A);
    LUA.inv(A);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
bool Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::retrieve
(
    const Foam::scalarField& phiq,
    scalarField& Rphiq
)
{
    bool retrieved(false);
    chemPointISAT<CompType, ThermoType>* phi0;

    // If the tree is not empty
    if (chemisTree_.size())
    {
        chemisTree_.binaryTreeSearch(phiq, chemisTree_.root(), phi0);

        // lastSearch keeps track of the chemPoint we obtain by the regular
        // binary tree search
        lastSearch_ = phi0;
        if (phi0->inEOA(phiq))
        {
            retrieved = true;
        }
        // After a successful secondarySearch, phi0 store a pointer to the
        // found chemPoint
        else if (chemisTree_.secondaryBTSearch(phiq, phi0))
        {
            retrieved = true;
        }
        else if (MRURetrieve_)
        {
            typename SLList
            <
                chemPointISAT<CompType, ThermoType>*
            >::iterator iter = MRUList_.begin();

            for ( ; iter != MRUList_.end(); ++iter)
            {
                phi0 = iter();
                if (phi0->inEOA(phiq))
                {
                    retrieved = true;
                    break;
                }
            }
        }
    }
    // The tree is empty, retrieved is still false
    else
    {
        // There is no chempoints that we can try to grow
        lastSearch_ = nullptr;
    }

    if (retrieved)
    {
        phi0->increaseNumRetrieve();
        scalar elapsedTimeSteps =
            this->chemistry_.timeSteps() - phi0->timeTag();

        // Raise a flag when the chemPoint has been used more than the allowed
        // number of time steps
        if (elapsedTimeSteps > chPMaxLifeTime_ && !phi0->toRemove())
        {
            cleaningRequired_ = true;
            phi0->toRemove() = true;
        }
        lastSearch_->lastTimeUsed() = this->chemistry_.timeSteps();
        addToMRU(phi0);
        calcNewC(phi0,phiq, Rphiq);
        nRetrieved_++;
        return true;
    }
    else
    {
        // This point is reached when every retrieve trials have failed
        // or if the tree is empty
        return false;
    }
}


template<class CompType, class ThermoType>
Foam::label Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::add
(
    const scalarField& phiq,
    const scalarField& Rphiq,
    const scalar rho,
    const scalar deltaT
)
{
    label growthOrAddFlag = 1;
    // If lastSearch_ holds a valid pointer to a chemPoint AND the growPoints_
    // option is on, the code first tries to grow the point hold by lastSearch_
    if (lastSearch_ && growPoints_)
    {
        if (grow(lastSearch_,phiq, Rphiq))
        {
            nGrowth_++;
            growthOrAddFlag = 0;
            // the structure of the tree is not modified, return false
            return growthOrAddFlag;
        }
    }

    // If the code reach this point, it is either because lastSearch_ is not
    // valid, OR because growPoints_ is not on, OR because the grow operation
    // has failed. In the three cases, a new point is added to the tree.
    if (chemisTree().isFull())
    {
        // If cleanAndBalance operation do not result in a reduction of the tree
        // size, the last possibility is to delete completely the tree.
        // It can be partially rebuild with the MRU list if this is used.
        if (!cleanAndBalance())
        {
            DynamicList<chemPointISAT<CompType, ThermoType>*> tempList;
            if (maxMRUSize_>0)
            {
                // Create a copy of each chemPointISAT of the MRUList_ before
                // they are deleted
                typename SLList
                <
                    chemPointISAT<CompType, ThermoType>*
                >::iterator iter = MRUList_.begin();
                for ( ; iter != MRUList_.end(); ++iter)
                {
                    tempList.append
                    (
                        new chemPointISAT<CompType, ThermoType>(*iter())
                    );
                }
            }
            chemisTree().clear();

            // Pointers to chemPoint are not valid anymore, clear the list
            MRUList_.clear();

            // Construct the tree without giving a reference to attach to it
            // since the structure has been completely discarded
            chemPointISAT<CompType, ThermoType>* nulPhi = 0;
            forAll(tempList, i)
            {
                chemisTree().insertNewLeaf
                (
                     tempList[i]->phi(),
                     tempList[i]->Rphi(),
                     tempList[i]->A(),
                     scaleFactor(),
                     this->tolerance(),
                     scaleFactor_.size(),
                     nulPhi
                );
                deleteDemandDrivenData(tempList[i]);
            }
        }

        // The structure has been changed, it will force the binary tree to
        // perform a new search and find the most appropriate point still stored
        lastSearch_ = nullptr;
    }

    // Compute the A matrix needed to store the chemPoint.
    label ASize = this->chemistry_.nEqns() + nAdditionalEqns_ - 2;
    scalarSquareMatrix A(ASize, Zero);
    computeA(A, Rphiq, rho, deltaT);

    chemisTree().insertNewLeaf
    (
        phiq,
        Rphiq,
        A,
        scaleFactor(),
        this->tolerance(),
        scaleFactor_.size(),
        lastSearch_ // lastSearch_ may be nullptr (handled by binaryTree)
    );

    nAdd_++;

    return growthOrAddFlag;
}


template<class CompType, class ThermoType>
void
Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::writePerformance()
{
    if (this->log())
    {
        nRetrievedFile_()
            << runTime_.timeOutputValue() << "    " << nRetrieved_ << endl;
        nRetrieved_ = 0;

        nGrowthFile_()
            << runTime_.timeOutputValue() << "    " << nGrowth_ << endl;
        nGrowth_ = 0;

        nAddFile_()
            << runTime_.timeOutputValue() << "    " << nAdd_ << endl;
        nAdd_ = 0;

        sizeFile_()
            << runTime_.timeOutputValue() << "    " << this->size() << endl;
    }
}


// ************************************************************************* //
