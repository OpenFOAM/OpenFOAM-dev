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

#include "binaryTree.H"
#include "SortableList.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

bool Foam::binaryTree::inSubTree
(
    const scalarField& phiq,
    binaryNode* y,
    chemPointISAT* x
)
{
    if ((n2ndSearch_ < max2ndSearch_) && (y!=nullptr))
    {
        scalar vPhi = 0;
        const scalarField& v = y->v();
        const scalar a = y->a();
        // compute v*phi
        for (label i=0; i<phiq.size(); i++)
        {
            vPhi += phiq[i]*v[i];
        }
        if (vPhi <= a)// on the left side of the node
        {
            if (y->nodeLeft() == nullptr)// left is a chemPoint
            {
                n2ndSearch_++;
                if (y->leafLeft()->inEOA(phiq))
                {
                    x = y->leafLeft();
                    return true;
                }
            }
            else // the left side is a node
            {
                if (inSubTree(phiq, y->nodeLeft(),x))
                {
                    return true;
                }
            }

            // not on the left side, try the right side
            if ((n2ndSearch_ < max2ndSearch_) && y->nodeRight() == nullptr)
            {
                n2ndSearch_++;
                // we reach the end of the subTree we can return the result
                if (y->leafRight()->inEOA(phiq))
                {
                    x = y->leafRight();
                    return true;
                }
                else
                {
                    x = nullptr;
                    return false;
                }
            }
            else // test for n2ndSearch is done in the call of inSubTree
            {
                return inSubTree(phiq, y->nodeRight(),x);
            }
        }
        else // on right side (symmetric of above)
        {
            if (y->nodeRight() == nullptr)
            {
                n2ndSearch_++;
                if (y->leafRight()->inEOA(phiq))
                {
                    return true;
                }
            }
            else // the right side is a node
            {
                if (inSubTree(phiq, y->nodeRight(),x))
                {
                    x = y->leafRight();
                    return true;
                }
            }
            // if we reach this point, the retrieve has
            // failed on the right side, explore the left side
            if ((n2ndSearch_ < max2ndSearch_) && y->nodeLeft() == nullptr)
            {
                n2ndSearch_++;
                if (y->leafLeft()->inEOA(phiq))
                {
                    x = y->leafLeft();
                    return true;
                }
                else
                {
                    x = nullptr;
                    return false;
                }
            }
            else
            {
                return inSubTree(phiq, y->nodeLeft(),x);
            }
        }
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binaryTree::binaryTree
(
    chemistryTabulationMethods::ISAT& table,
    dictionary coeffsDict
)
:
    table_(table),
    root_(nullptr),
    maxNLeafs_(coeffsDict.lookup<label>("maxNLeafs")),
    size_(0),
    n2ndSearch_(0),
    max2ndSearch_(coeffsDict.lookupOrDefault("max2ndSearch",0)),
    coeffsDict_(coeffsDict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::binaryTree::insertNewLeaf
(
    const scalarField& phiq,
    const scalarField& Rphiq,
    const scalarSquareMatrix& A,
    const scalarField& scaleFactor,
    const scalar& epsTol,
    const label nCols,
    const label nActive,
    chemPointISAT*& phi0
)
{
    if (size_ == 0) // no points are stored
    {
        // create an empty binary node and point root_ to it
        root_ = new binaryNode();
        // create the new chemPoint which holds the composition point
        // phiq and the data to initialise the EOA
        chemPointISAT* newChemPoint =
            new chemPointISAT
            (
                table_,
                phiq,
                Rphiq,
                A,
                scaleFactor,
                epsTol,
                nCols,
                nActive,
                coeffsDict_,
                root_
            );
        root_->leafLeft() = newChemPoint;
    }
    else // at least one point stored
    {
        // no reference chemPoint, a BT search is required
        if (phi0 == nullptr)
        {
            binaryTreeSearch(phiq, root_,phi0);
        }
        // access to the parent node of the chemPoint
        binaryNode* parentNode = phi0->node();

        // create the new chemPoint which holds the composition point
        // phiq and the data to initialise the EOA
        chemPointISAT* newChemPoint =
            new chemPointISAT
            (
                table_,
                phiq,
                Rphiq,
                A,
                scaleFactor,
                epsTol,
                nCols,
                nActive,
                coeffsDict_
            );
        // insert new node on the parent node in the position of the
        // previously stored leaf (phi0)
        // the new node contains phi0 on the left and phiq on the right
        // the hyper plane is computed in the binaryNode constructor
        binaryNode* newNode;
        if (size_>1)
        {
            newNode = new binaryNode(phi0, newChemPoint, parentNode);
            // make the parent of phi0 point to the newly created node
            insertNode(phi0, newNode);
        }
        else // size_ == 1 (because not equal to 0)
        {
            // when size is 1, the binaryNode is without hyperplane
            deleteDemandDrivenData(root_);
            newNode = new binaryNode(phi0, newChemPoint, nullptr);
            root_ = newNode;
        }

        phi0->node() = newNode;
        newChemPoint->node()=newNode;
    }
    size_++;
}


bool Foam::binaryTree::secondaryBTSearch
(
    const scalarField& phiq,
    chemPointISAT*& x
)
{
    // initialise n2ndSearch_
    n2ndSearch_ = 0;
    if ((n2ndSearch_ < max2ndSearch_) && (size_ > 1))
    {
        chemPointISAT* xS = chemPSibling(x);
        if (xS != nullptr)
        {
            n2ndSearch_++;
            if (xS->inEOA(phiq))
            {
                x = xS;
                return true;
            }
        }
        else if (inSubTree(phiq, nodeSibling(x),x))
        {
            return true;
        }
        // if we reach this point, no leafs were found at this depth or lower
        // we move upward in the tree
        binaryNode* y = x->node();
        while((y->parent()!= nullptr) && (n2ndSearch_ < max2ndSearch_))
        {
            xS = chemPSibling(y);
            if (xS != nullptr)
            {
                n2ndSearch_++;
                if (xS->inEOA(phiq))
                {
                    x=xS;
                    return true;
                }
            }
            else if (inSubTree(phiq, nodeSibling(y),x))
            {
                return true;
            }
            y = y->parent();
        }
        // if we reach this point it is either because
        // we did not find another covering EOA in the entire tree or
        // we reach the maximum number of secondary search
        return false;
    }
    else
    {
        return false;
    }
}


void Foam::binaryTree::deleteLeaf(chemPointISAT*& phi0)
{
    if (size_ == 1) // only one point is stored
    {
        deleteDemandDrivenData(phi0);
        deleteDemandDrivenData(root_);
    }
    else if (size_ > 1)
    {
        binaryNode* z = phi0->node();
        binaryNode* x;
        chemPointISAT* siblingPhi0 = chemPSibling(phi0);

        if (siblingPhi0 != nullptr)// the sibling of phi0 is a chemPoint
        {
            // z was root (only two chemPoints in the tree)
            if (z->parent() == nullptr)
            {
                root_ = new binaryNode();
                root_->leafLeft()=siblingPhi0;
                siblingPhi0->node()=root_;
            }
            else if (z == z->parent()->nodeLeft())
            {
                z->parent()->leafLeft() = siblingPhi0;
                z->parent()->nodeLeft() = nullptr;
                siblingPhi0->node() = z->parent();
            }
            else if (z == z->parent()->nodeRight())
            {
                z->parent()->leafRight() = siblingPhi0;
                z->parent()->nodeRight() = nullptr;
                siblingPhi0->node() = z->parent();
            }
            else
            {
                FatalErrorInFunction
                    << "wrong addressing of the initial leaf"
                    << exit(FatalError);
            }
        }
        else
        {
            x = nodeSibling(phi0);
            if (x !=nullptr)
            {
                transplant(z, x);
            }
            else
            {
                FatalErrorInFunction
                    << "inconsistent structure of the tree, no leaf and no node"
                    << exit(FatalError);
            }
        }
        deleteDemandDrivenData(phi0);
        deleteDemandDrivenData(z);
    }
    size_--;
}


void Foam::binaryTree::balance()
{
    //1) walk through the entire tree by starting with the tree's most left
    // chemPoint
    chemPointISAT* x = treeMin();
    List<chemPointISAT*> chemPoints(size_);
    label chemPointISATi=0;
    while (x != nullptr)
    {
        chemPoints[chemPointISATi++] = x;
        x = treeSuccessor(x);
    }

    //2) compute the mean composition
    scalarField mean(treeMin()->phi().size(), Zero);
    forAll(chemPoints, j)
    {
        const scalarField& phij = chemPoints[j]->phi();
        mean += phij;
    }
    mean /= size_;

    //3) compute the variance for each space direction
    List<scalar> variance(treeMin()->phi().size(), Zero);
    forAll(chemPoints, j)
    {
        const scalarField& phij = chemPoints[j]->phi();
        forAll(variance, vi)
        {
            variance[vi] += sqr(phij[vi]-mean[vi]);
        }
    }

    //4) analyze what is the direction of the maximal variance
    scalar maxVariance(-1.0);
    label maxDir(-1);
    forAll(variance, vi)
    {
        if (maxVariance < variance[vi])
        {
            maxVariance = variance[vi];
            maxDir = vi;
        }
    }

    // maxDir indicates the direction of maximum variance
    // we create the new root node by taking the two extreme points
    // in this direction if these extreme points were not deleted in the
    // cleaning that come before the balance function they are still important
    // and the tree should therefore take them into account
    SortableList<scalar> phiMaxDir(chemPoints.size(),0.0);
    forAll(chemPoints, j)
    {
        phiMaxDir[j] = chemPoints[j]->phi()[maxDir];
    }

    phiMaxDir.sort();
    // delete reference to all node since the tree is reshaped
    deleteAllNode();
    root_ = nullptr;

    // add the node for the two extremum
    binaryNode* newNode = new binaryNode
    (
        chemPoints[phiMaxDir.indices()[0]],
        chemPoints[phiMaxDir.indices()[phiMaxDir.size()-1]],
        nullptr
    );
    root_ = newNode;

    chemPoints[phiMaxDir.indices()[0]]->node() = newNode;
    chemPoints[phiMaxDir.indices()[phiMaxDir.size()-1]]->node() = newNode;

    for (label cpi=1; cpi<chemPoints.size()-1; cpi++)
    {
        chemPointISAT* phi0;
        binaryTreeSearch
        (
            chemPoints[phiMaxDir.indices()[cpi]]->phi(),
            root_,
            phi0
        );
        // add the chemPoint
        binaryNode* nodeToAdd = new binaryNode
        (
            phi0,
            chemPoints[phiMaxDir.indices()[cpi]],
            phi0->node()
        );

        // make the parent of phi0 point to the newly created node
        insertNode(phi0, nodeToAdd);
        phi0->node() = nodeToAdd;
        chemPoints[phiMaxDir.indices()[cpi]]->node() = nodeToAdd;
    }
}


Foam::chemPointISAT* Foam::binaryTree::treeSuccessor(chemPointISAT* x)
{
    if (size_>1)
    {
        if (x == x->node()->leafLeft())
        {
            if (x->node()->nodeRight() == nullptr)
            {
                return x->node()->leafRight();
            }
            else
            {
                return treeMin(x->node()->nodeRight());
            }
        }
        else if (x == x->node()->leafRight())
        {
            binaryNode* y = x->node();
            while((y->parent() !=nullptr))
            {
                if (y == y->parent()->nodeLeft())
                {
                    if (y->parent()->nodeRight() == nullptr)
                    {
                        return y->parent()->leafRight();
                    }
                    else
                    {
                        return treeMin(y->parent()->nodeRight());
                    }
                }
                y=y->parent();
            }
            // when we reach this point, y points to the root and
            // never entered in the if loop (coming from the right)
            // so we are at the tree maximum and there is no successor
            return nullptr;
        }
        else
        {
            FatalErrorInFunction
                << "inconsistent structure of the tree, no leaf and no node"
                << exit(FatalError);
            return nullptr;
        }
    }

    return nullptr;
}


// ************************************************************************* //
