/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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
#include "odeChemistryModel.H"
#include "LUscalarMatrix.H"
#include "addToRunTimeSelectionTable.H"


/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace chemistryTabulationMethods
{
    defineTypeNameAndDebug(ISAT, 0);
    addToRunTimeSelectionTable(chemistryTabulationMethod, ISAT, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chemistryTabulationMethods::ISAT::ISAT
(
    const dictionary& chemistryProperties,
    const odeChemistryModel& chemistry
)
:
    chemistryTabulationMethod
    (
        chemistryProperties,
        chemistry
    ),
    coeffsDict_(chemistryProperties.subDict("tabulation")),
    chemistry_(chemistry),
    log_(coeffsDict_.lookupOrDefault<Switch>("log", false)),
    reduction_(chemistry_.reduction()),
    chemisTree_(*this, coeffsDict_),
    scaleFactor_(chemistry.nEqns() + 1, 1),
    runTime_(chemistry.time()),
    timeSteps_(0),
    chPMaxLifeTime_
    (
        coeffsDict_.lookupOrDefault("chPMaxLifeTime", INT_MAX)
    ),
    maxGrowth_(coeffsDict_.lookupOrDefault("maxGrowth", INT_MAX)),
    checkEntireTreeInterval_
    (
        coeffsDict_.lookupOrDefault("checkEntireTreeInterval", INT_MAX)
    ),
    maxDepthFactor_
    (
        coeffsDict_.lookupOrDefault
        (
            "maxDepthFactor",
            (chemisTree_.maxNLeafs() - 1)
           /(Foam::log(scalar(chemisTree_.maxNLeafs()))/Foam::log(2.0))
        )
    ),
    minBalanceThreshold_
    (
        coeffsDict_.lookupOrDefault
        (
            "minBalanceThreshold", 0.1*chemisTree_.maxNLeafs()
        )
    ),
    MRURetrieve_(coeffsDict_.lookupOrDefault("MRURetrieve", false)),
    maxMRUSize_(coeffsDict_.lookupOrDefault("maxMRUSize", 0)),
    lastSearch_(nullptr),
    growPoints_(coeffsDict_.lookupOrDefault("growPoints", true)),
    tolerance_(coeffsDict_.lookupOrDefault("tolerance", 1e-4)),
    nRetrieved_(0),
    nGrowth_(0),
    nAdd_(0),
    addNewLeafCpuTime_(0),
    growCpuTime_(0),
    searchISATCpuTime_(0),
    tabulationResults_
    (
        IOobject
        (
            chemistry.thermo().phasePropertyName("TabulationResults"),
            chemistry.time().name(),
            chemistry.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        chemistry.mesh(),
        scalar(0)
    ),

    cleaningRequired_(false)
{
    dictionary scaleDict(coeffsDict_.subDict("scaleFactor"));
    label Ysize = chemistry_.Y().size();
    scalar otherScaleFactor = scaleDict.lookup<scalar>("otherSpecies");
    for (label i=0; i<Ysize; i++)
    {
        if (!scaleDict.found(chemistry_.Y()[i].member()))
        {
            scaleFactor_[i] = otherScaleFactor;
        }
        else
        {
            scaleFactor_[i] =
                scaleDict.lookup<scalar>(chemistry_.Y()[i].member());
        }
    }
    scaleFactor_[Ysize] = scaleDict.lookup<scalar>("Temperature");
    scaleFactor_[Ysize + 1] = scaleDict.lookup<scalar>("Pressure");
    scaleFactor_[Ysize + 2] = scaleDict.lookup<scalar>("deltaT");

    if (log_)
    {
        nRetrievedFile_ = chemistry.logFile("found_isat.out");
        nGrowthFile_ = chemistry.logFile("growth_isat.out");
        nAddFile_ = chemistry.logFile("add_isat.out");
        sizeFile_ = chemistry.logFile("size_isat.out");

        cpuAddFile_ = chemistry.logFile("cpu_add.out");
        cpuGrowFile_ = chemistry.logFile("cpu_grow.out");
        cpuRetrieveFile_ = chemistry.logFile("cpu_retrieve.out");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::chemistryTabulationMethods::ISAT::~ISAT()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::chemistryTabulationMethods::ISAT::addToMRU
(
    chemPointISAT* phi0
)
{
    if (maxMRUSize_ > 0 && MRURetrieve_)
    {
        // First search if the chemPoint is already in the list
        bool isInList = false;
        typename SLList <chemPointISAT*>::iterator iter =
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
                MRUList_.remove(iter);
                MRUList_.insert(phi0);
            }
            else
            {
                MRUList_.insert(phi0);
            }
        }
    }
}


void Foam::chemistryTabulationMethods::ISAT::calcNewC
(
    chemPointISAT* phi0,
    const scalarField& phiq,
    scalarField& Rphiq
)
{
    const label nEqns = chemistry_.nEqns(); // Species, T, p
    const List<label>& completeToSimplified = phi0->completeToSimplifiedIndex();

    const scalarField dphi(phiq - phi0->phi());
    const scalarSquareMatrix& gradientsMatrix = phi0->A();

    // Linear extrapolation:
    //
    //     Rphiq[i] = Rphi0[i] + (dR_i/dphi_j)*dphi[j]
    //              = Rphi0[i] + A(i, j)*dphi[j]
    //

    Rphiq = phi0->Rphi();
    for (label i=0; i<nEqns - 2; i++)
    {
        if (reduction_)
        {
            const label si = completeToSimplified[i];

            if (si != -1)
            {
                // If specie is active, or T or p, then extrapolate using the
                // gradients matrix
                for (label j=0; j<nEqns + 1; j++)
                {
                    const label sj =
                        j < nEqns - 2
                      ? completeToSimplified[j]
                      : j - (nEqns - 2) + phi0->nActive();

                    if (sj != -1)
                    {
                        Rphiq[i] += gradientsMatrix(si, sj)*dphi[j];
                    }
                }
            }
            else
            {
                // If specie is inactive then use the tabulated value directly
                Rphiq[i] += dphi[i];
            }
        }
        else
        {
            // Extrapolate using the gradients matrix
            for (label j=0; j<nEqns + 1; j++)
            {
                Rphiq[i] += gradientsMatrix(i, j)*dphi[j];
            }
        }

        // Clip
        Rphiq[i] = max(0, Rphiq[i]);
    }
}


bool Foam::chemistryTabulationMethods::ISAT::grow
(
    chemPointISAT* phi0,
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


bool Foam::chemistryTabulationMethods::ISAT::cleanAndBalance()
{
    bool treeModified(false);

    // Check all chemPoints to see if we need to delete some of the chemPoints
    // according to the elapsed time and number of growths
    chemPointISAT* x = chemisTree_.treeMin();
    while(x != nullptr)
    {
        chemPointISAT* xtmp = chemisTree_.treeSuccessor(x);

        const scalar elapsedTimeSteps = timeSteps() - x->timeTag();

        if ((elapsedTimeSteps > chPMaxLifeTime_) || (x->nGrowth() > maxGrowth_))
        {
            chemisTree_.deleteLeaf(x);
            treeModified = true;
        }
        x = xtmp;
    }

    MRUList_.clear();

    // Check if the tree should be balanced according to criterion:
    //  -the depth of the tree bigger than a*log2(size), log2(size) being the
    //      ideal depth (e.g. 4 leafs can be stored in a tree of depth 2)
    if
    (
        chemisTree_.size() > minBalanceThreshold_
     && chemisTree_.depth() >
        maxDepthFactor_*Foam::log(scalar(chemisTree_.size()))/Foam::log(2.0)
    )
    {
        chemisTree_.balance();
        treeModified = true;
    }

    // Return a bool to specify if the tree structure has been modified and is
    // now below the user specified limit (true if not full)
    return (treeModified && !chemisTree_.isFull());
}


void Foam::chemistryTabulationMethods::ISAT::computeA
(
    scalarSquareMatrix& A,
    const scalarField& Rphiq,
    const label li,
    const scalar dt
)
{
    const label nSpecie = chemistry_.nSpecie();

    scalarField Rphiqs(chemistry_.nEqns() + 1);
    for (label i=0; i<nSpecie; i++)
    {
        const label si = chemistry_.sToc(i);
        Rphiqs[i] = Rphiq[si];
    }
    Rphiqs[nSpecie] = Rphiq[Rphiq.size() - 3];
    Rphiqs[nSpecie + 1] = Rphiq[Rphiq.size() - 2];
    Rphiqs[nSpecie + 2] = Rphiq[Rphiq.size() - 1];

    // Aaa is computed implicitly,
    // A is given by A = C(psi0, t0+dt), where C is obtained through solving
    // d/dt C(psi0, t) = J(psi(t))C(psi0, t)
    // If we solve it implicitly:
    // (C(psi0, t0+dt) - C(psi0, t0))/dt = J(psi(t0+dt))C(psi0, t0+dt)
    // The Jacobian is thus computed according to the mapping
    // C(psi0,t0+dt)*(I-dt*J(psi(t0+dt))) = C(psi0, t0)
    // A = C(psi0,t0)/(I-dt*J(psi(t0+dt)))
    // where C(psi0,t0) = I
    scalarField dYTpdt(nSpecie + 2, Zero);
    chemistry_.jacobian(runTime_.value(), Rphiqs, li, dYTpdt, A);

    // Inverse of I - dt*J(psi(t0 + dt))
    for (label i=0; i<nSpecie + 2; i++)
    {
        for (label j=0; j<nSpecie + 2; j++)
        {
            A(i, j) *= -dt;
        }
        A(i, i) += 1;
    }
    A(nSpecie + 2, nSpecie + 2) = 1;
    LUscalarMatrix LUA(A);
    LUA.inv(A);

    // After inversion, lines of p and T are set to 0 except diagonal.  This
    // avoid skewness of the ellipsoid of accuracy and potential issues in the
    // binary tree.
    for (label i=0; i<nSpecie; i++)
    {
        A(nSpecie, i) = 0;
        A(nSpecie + 1, i) = 0;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::chemistryTabulationMethods::ISAT::retrieve
(
    const Foam::scalarField& phiq,
    scalarField& Rphiq
)
{
    if (log_)
    {
        cpuTime_.cpuTimeIncrement();
    }

    bool retrieved(false);
    chemPointISAT* phi0;

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
                chemPointISAT*
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
        const scalar elapsedTimeSteps = timeSteps() - phi0->timeTag();

        // Raise a flag when the chemPoint has been used more than the allowed
        // number of time steps
        if (elapsedTimeSteps > chPMaxLifeTime_ && !phi0->toRemove())
        {
            cleaningRequired_ = true;
            phi0->toRemove() = true;
        }
        lastSearch_->lastTimeUsed() = timeSteps();
        addToMRU(phi0);
        calcNewC(phi0, phiq, Rphiq);
        nRetrieved_++;
    }

    if (log_)
    {
        searchISATCpuTime_ += cpuTime_.cpuTimeIncrement();
    }

    return retrieved;
}


Foam::label Foam::chemistryTabulationMethods::ISAT::add
(
    const scalarField& phiq,
    const scalarField& Rphiq,
    const label nActive,
    const label li,
    const scalar deltaT
)
{
    if (log_)
    {
        cpuTime_.cpuTimeIncrement();
    }

    label growthOrAddFlag = 1;

    // If lastSearch_ holds a valid pointer to a chemPoint AND the growPoints_
    // option is on, the code first tries to grow the point hold by lastSearch_
    if (lastSearch_ && growPoints_)
    {
        if (grow(lastSearch_, phiq, Rphiq))
        {
            nGrowth_++;
            growthOrAddFlag = 0;
            addToMRU(lastSearch_);

            tabulationResults_[li] = 1;

            if (log_)
            {
                growCpuTime_ += cpuTime_.cpuTimeIncrement();
            }

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
            DynamicList<chemPointISAT*> tempList;
            if (maxMRUSize_>0)
            {
                // Create a copy of each chemPointISAT of the MRUList_ before
                // they are deleted
                typename SLList
                <
                    chemPointISAT*
                >::iterator iter = MRUList_.begin();
                for ( ; iter != MRUList_.end(); ++iter)
                {
                    tempList.append
                    (
                        new chemPointISAT(*iter())
                    );
                }
            }
            chemisTree().clear();

            // Pointers to chemPoint are not valid anymore, clear the list
            MRUList_.clear();

            // Construct the tree without giving a reference to attach to it
            // since the structure has been completely discarded
            chemPointISAT* nulPhi = 0;
            forAll(tempList, i)
            {
                chemisTree().insertNewLeaf
                (
                     tempList[i]->phi(),
                     tempList[i]->Rphi(),
                     tempList[i]->A(),
                     scaleFactor(),
                     tolerance_,
                     scaleFactor_.size(),
                     nActive,
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
    const label ASize = chemistry_.nEqns() + 1;
    scalarSquareMatrix A(ASize, Zero);
    computeA(A, Rphiq, li, deltaT);

    chemisTree().insertNewLeaf
    (
        phiq,
        Rphiq,
        A,
        scaleFactor(),
        tolerance_,
        scaleFactor_.size(),
        nActive,
        lastSearch_ // lastSearch_ may be nullptr (handled by binaryTree)
    );
    if (lastSearch_ != nullptr)
    {
        addToMRU(lastSearch_);
    }
    nAdd_++;

    tabulationResults_[li] = 0;

    if (log_)
    {
        addNewLeafCpuTime_ += cpuTime_.cpuTimeIncrement();
    }

    return growthOrAddFlag;
}


void Foam::chemistryTabulationMethods::ISAT::writePerformance()
{
    if (log_)
    {
        nRetrievedFile_()
            << runTime_.userTimeValue() << "    " << nRetrieved_ << endl;
        nRetrieved_ = 0;

        nGrowthFile_()
            << runTime_.userTimeValue() << "    " << nGrowth_ << endl;
        nGrowth_ = 0;

        nAddFile_()
            << runTime_.userTimeValue() << "    " << nAdd_ << endl;
        nAdd_ = 0;

        sizeFile_()
            << runTime_.userTimeValue() << "    "
            << chemisTree_.size() << endl;

        cpuRetrieveFile_()
            << runTime_.userTimeValue()
            << "    " << searchISATCpuTime_ << endl;
        searchISATCpuTime_ = 0;

        cpuGrowFile_()
            << runTime_.userTimeValue()
            << "    " << growCpuTime_ << endl;
        growCpuTime_ = 0;

        cpuAddFile_()
            << runTime_.userTimeValue()
            << "    " << addNewLeafCpuTime_ << endl;
        addNewLeafCpuTime_ = 0;
    }
}


void Foam::chemistryTabulationMethods::ISAT::reset()
{
    // Increment counter of time-step
    timeSteps_++;

    forAll(tabulationResults_, i)
    {
        tabulationResults_[i] = 2;
    }
}


bool Foam::chemistryTabulationMethods::ISAT::update()
{
    bool updated = cleanAndBalance();
    writePerformance();
    return updated;
}


// ************************************************************************* //
