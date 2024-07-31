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

Class
    Foam::binaryTree

Description

    Data storage of the chemistryOnLineLibrary according to a binary
    tree structure.

            0 (root node)
         /     \
        0       0
      /   \   /   \
     L     R L     0
    / \
   L   R

    L: leafLeft_
    R: leafRight_

\*---------------------------------------------------------------------------*/

#ifndef binaryTree_H
#define binaryTree_H

#include "binaryNode.H"
#include "chemPointISAT.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class binaryTree Declaration
\*---------------------------------------------------------------------------*/

class binaryTree
{
    // Private Data

        //- Reference to the ISAT table
        chemistryTabulationMethods::ISAT& table_;

        //- Root node of the binary tree
        binaryNode *root_;

        //- Maximum number of elements in the binary tree
        label maxNLeafs_;

        //- Size of the BST (= number of chemPoint stored)
        label size_;

        //- Secondary retrieve search variables
        label n2ndSearch_;

        //- Secondary retrieve search variables
        label max2ndSearch_;

        dictionary coeffDict_;


    // Private Member Functions

        //- Insert new node at the position of phi0. phi0 should be already
        //  attached to another node or the pointer to it will be lost.
        inline void insertNode(chemPointISAT*& phi0, binaryNode*& newNode);

        //- Perform a search in the subtree starting from the subtree node y.
        //  This search continues to use the hyperplane to walk the tree.
        //  If covering EOA is found return true and x points to the chemPoint.
        bool inSubTree
        (
            const scalarField& phiq,
            binaryNode* y,
            chemPointISAT* x
        );

        inline void deleteSubTree(binaryNode* subTreeRoot);

        inline void deleteSubTree();

        //- Replace the binaryNode u with v
        inline void transplant(binaryNode* u, binaryNode* v);

        inline chemPointISAT* chemPSibling(binaryNode* y);

        inline chemPointISAT* chemPSibling(chemPointISAT* x);

        inline binaryNode* nodeSibling(binaryNode* y);

        inline binaryNode* nodeSibling(chemPointISAT* x);

        inline void deleteAllNode(binaryNode* subTreeRoot);


public:

    //- Constructors

        //- Construct from dictionary and chemistryOnLineLibrary
        binaryTree
        (
            chemistryTabulationMethods::ISAT& table,
            dictionary coeffDict
        );


    //- Member functions

        inline label size() const;

        //- Computes iteratively the depth of the subTree
        inline label depth(binaryNode* subTreeRoot) const;

        inline label depth() const;

        inline binaryNode* root();

        inline label maxNLeafs() const;

        // Insert a new leaf starting from the parent node of phi0
        // Parameters: phi0 the leaf to replace by a node
        // phiq the new composition to store
        // Rphiq the mapping of the new composition point
        // A the mapping gradient matrix
        // B the matrix used to initialise the EOA
        // nCols the size of the matrix
        // Returns: void
        // Description :
        //1) Create a new leaf with the data to initialise the EOA and to
        // retrieve the mapping by linear interpolation (the EOA is
        // initialise in the chemPoint constructor)
        //2) Get the parent node of phi0 and connect a new node in place of the
        // leaf of phi0. This new node is constructed with phi0 on the left
        // and phiq on the right (the hyperplane is computed inside the
        // binaryNode constructor)
        void insertNewLeaf
        (
            const scalarField& phiq,
            const scalarField& Rphiq,
            const scalarSquareMatrix& A,
            const scalarField& scaleFactor,
            const scalar& epsTol,
            const label nCols,
            const label nActive,
            chemPointISAT*& phi0
        );

        // Search the binaryTree until the nearest leaf of a specified
        // leaf is found.
        void binaryTreeSearch
        (
            const scalarField& phiq,
            binaryNode* node,
            chemPointISAT*& nearest
        );

        // Perform a secondary binary tree search starting from a failed
        // chemPoint x, with a depth-first search algorithm
        // If another candidate is found return true and x points to the chemP
        bool secondaryBTSearch(const scalarField& phiq, chemPointISAT*& x);

        //- Delete a leaf from the binary tree and reshape the binary tree for
        //  the following binary tree search
        //  Return the index in the nodeList of the removed node
        //  (-1 when no node)
        void deleteLeaf(chemPointISAT*& phi0);

        //- Cheap balance function
        //  This function just roughly separate the space in two parts
        //  with a hyperplane which separate the two extreme chemPoint in the
        //  direction of the maximum the variance
        //  Then, it repopulate the tree with this hyperplane stored at the root
        //  and by inserting the chemPoint in increasing order of value in that
        //  direction
        void balance();

        inline void deleteAllNode()
        {
            deleteAllNode(root_);
        }

        inline chemPointISAT* treeMin(binaryNode* subTreeRoot);

        inline chemPointISAT* treeMin()
        {
            return treeMin(root_);
        }

        chemPointISAT* treeSuccessor(chemPointISAT* x);

        //- Removes every entries of the tree and delete the associated objects
        inline void clear();

        //- ListFull
        inline bool isFull();

        inline void resetNumRetrieve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "binaryTreeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
