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

#include "demandDrivenData.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

inline void Foam::binaryTree::insertNode
(
    chemPointISAT*& phi0,
    binaryNode*& newNode
)
{
    if (phi0 == phi0->node()->leafRight())// phi0 is on the right
    {
        phi0->node()->leafRight() = nullptr;
        phi0->node()->nodeRight() = newNode;
        return;
    }
    else if (phi0 == phi0->node()->leafLeft())// phi0 is on the left
    {
        phi0->node()->leafLeft() = nullptr;
        phi0->node()->nodeLeft() = newNode;
        return;

    }

    // if we reach this point, there is an issue with the addressing
    FatalErrorInFunction
        << "trying to insert a node with a wrong pointer to a chemPoint"
        << exit(FatalError);
}


inline void Foam::binaryTree::deleteSubTree(binaryNode* subTreeRoot)
{
    if (subTreeRoot != nullptr)
    {
        deleteDemandDrivenData(subTreeRoot->leafLeft());
        deleteDemandDrivenData(subTreeRoot->leafRight());
        deleteSubTree(subTreeRoot->nodeLeft());
        deleteSubTree(subTreeRoot->nodeRight());
        deleteDemandDrivenData(subTreeRoot);
    }
}


inline void Foam::binaryTree::deleteSubTree()
{
    deleteSubTree(root_);
}


inline void Foam::binaryTree::transplant(binaryNode* u, binaryNode* v)
{
    if (v != nullptr)
    {
        // u is root_
        if (u->parent() == nullptr)
        {
            root_ = v;
        }
        // u is on the left of its parent
        else if (u == u->parent()->nodeLeft())
        {
            u->parent()->nodeLeft() = v;
        }
        // u is ont the right of its parent
        else if (u == u->parent()->nodeRight())
        {
            u->parent()->nodeRight() = v;
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial node"
                << exit(FatalError);
        }
        v->parent() = u->parent();
    }
    else
    {
        FatalErrorInFunction
            << "trying to transplant a nullptr node"
            << exit(FatalError);
    }
}


inline Foam::chemPointISAT* Foam::binaryTree::chemPSibling(binaryNode* y)
{
    if (y->parent() != nullptr)
    {
        if (y == y->parent()->nodeLeft())// y is on the left, return right side
        {
            // might return nullptr if the right side is a node
            return y->parent()->leafRight();
        }
        else if (y == y->parent()->nodeRight())
        {
            return y->parent()->leafLeft();
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial node"
                << exit(FatalError);
            return nullptr;
        }
    }

    // the binaryNode is root_ and has no sibling
    return nullptr;
}


inline Foam::chemPointISAT* Foam::binaryTree::chemPSibling(chemPointISAT* x)
{
    if (size_>1)
    {
        if (x == x->node()->leafLeft())
        {
            // x is on the left, return right side
            // might return nullptr if the right side is a node
            return x->node()->leafRight();
        }
        else if (x == x->node()->leafRight())
        {
            // x is on the right, return left side
            return x->node()->leafLeft();
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial leaf"
                << exit(FatalError);
            return nullptr;
        }
    }

    // there is only one leaf attached to the root_, no sibling
    return nullptr;
}


inline Foam::binaryNode* Foam::binaryTree::nodeSibling(binaryNode* y)
{
    if (y->parent()!=nullptr)
    {
        if (y == y->parent()->nodeLeft())
        {
            // y is on the left, return right side
            return y->parent()->nodeRight();
        }
        else if (y == y->parent()->nodeRight())
        {
            return y->parent()->nodeLeft();
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial node"
                << exit(FatalError);
            return nullptr;
        }
    }
    return nullptr;
}


Foam::binaryNode* Foam::binaryTree::nodeSibling(chemPointISAT* x)
{
    if (size_>1)
    {
        if (x == x->node()->leafLeft())
        {
            // x is on the left, return right side
            return x->node()->nodeRight();
        }
        else if (x == x->node()->leafRight())
        {
            // x is on the right, return left side
            return x->node()->nodeLeft();
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial leaf"
                << exit(FatalError);
            return nullptr;
        }
    }
    return nullptr;
}


inline void Foam::binaryTree::deleteAllNode(binaryNode* subTreeRoot)
{
    if (subTreeRoot != nullptr)
    {
        deleteAllNode(subTreeRoot->nodeLeft());
        deleteAllNode(subTreeRoot->nodeRight());
        deleteDemandDrivenData(subTreeRoot);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::binaryTree::size() const
{
    return size_;
}


inline Foam::label Foam::binaryTree::depth(binaryNode* subTreeRoot) const
{
    // when we reach the leaf, we return 0
    if (subTreeRoot == nullptr)
    {
        return 0;
    }
    else
    {
        return
            1
          + max
            (
                depth(subTreeRoot->nodeLeft()),
                depth(subTreeRoot->nodeRight())
            );
    }
}


inline Foam::label Foam::binaryTree::depth() const
{
    return depth(root_);
}


inline Foam::binaryNode* Foam::binaryTree::root()
{
    return root_;
}


inline Foam::label Foam::binaryTree::maxNLeafs() const
{
    return maxNLeafs_;
}


inline void Foam::binaryTree::binaryTreeSearch
(
    const scalarField& phiq,
    binaryNode* node,
    chemPointISAT*& nearest
)
{
    if (size_ > 1)
    {
        scalar vPhi=0.0;
        const scalarField& v = node->v();
        const scalar& a = node->a();
        // compute v*phi
        for (label i=0; i<phiq.size(); i++) vPhi += phiq[i]*v[i];


        if (vPhi > a) // on right side (side of the newly added point)
        {
            if (node->nodeRight()!=nullptr)
            {
                binaryTreeSearch(phiq, node->nodeRight(), nearest);
            }
            else // the terminal node is reached, store leaf on the right
            {
                nearest = node->leafRight();
                return;
            }
        }
        else // on left side (side of the previously stored point)
        {
            if (node->nodeLeft()!=nullptr)
            {
                binaryTreeSearch(phiq, node->nodeLeft(), nearest);
            }
            else // the terminal node is reached, return element on right
            {
                nearest = node->leafLeft();
                return;
            }
        }
    }
    // only one point stored (left element of the root)
    else if (size_ == 1)
    {
        nearest = root_->leafLeft();
    }
    else // no point stored
    {
        nearest = nullptr;
    }
}


inline Foam::chemPointISAT* Foam::binaryTree::treeMin(binaryNode* subTreeRoot)
{
    if (subTreeRoot!=nullptr)
    {
        while(subTreeRoot->nodeLeft() != nullptr)
        {
            subTreeRoot = subTreeRoot->nodeLeft();
        }
        return subTreeRoot->leafLeft();
    }
    else
    {
        return nullptr;
    }
}


inline void Foam::binaryTree::clear()
{
    // Recursively delete the element in the subTree
    deleteSubTree();

    // Reset root node (should already be nullptr)
    root_ = nullptr;

    // Reset size_
    size_ = 0;
}


inline bool Foam::binaryTree::isFull()
{
    return size_ >= maxNLeafs_;
}


inline void Foam::binaryTree::resetNumRetrieve()
{
    // Has to go along each chemPointISAT of the tree
    if (size_ > 0)
    {
        // First finds the first leaf
        chemPointISAT* chP0 = treeMin();
        chP0->resetNumRetrieve();

        // Then go to each successor
        chemPointISAT* nextChP = treeSuccessor(chP0);
        while (nextChP != nullptr)
        {
            nextChP->resetNumRetrieve();
            nextChP = treeSuccessor(nextChP);
        }
    }
}


// ************************************************************************* //
