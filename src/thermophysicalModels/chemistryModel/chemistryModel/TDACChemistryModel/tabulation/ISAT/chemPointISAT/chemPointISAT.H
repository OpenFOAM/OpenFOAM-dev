/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

Class
    Foam::chemPointISAT

Description
    Leaf of the binary tree.
    The chemPoint stores the composition 'phi', the mapping of this
    composition Rphi, the mapping gradient matrix A and the matrix describing
    the Ellipsoid Of Accuracy (EOA).

  1)When the chemPoint is created the region of accuracy is approximated by
    an ellipsoid E centered in 'phi' (obtained with the constant):
        E = {x| ||L^T.(x-phi)|| <= 1},
    with x a point in the composition space and L^T the transpose of an upper
    triangular matrix describing the EOA (see below: "Computation of L" ).

  2)To RETRIEVE the mapping from the chemPoint phi, the query point phiq has to
    be in the EOA of phi. It follows that, dphi=phiq-phi and to test if phiq
    is in the ellipsoid there are two methods. First, compare r=||dphi|| with
    rmin and rmax. If r < rmin, phiq is in the EOA. If r > rmax, phiq is out of
    the EOA. This operations is O(completeSpaceSize) and is performed first.
    If rmin < r < rmax, then the second method is used:
        ||L^T.dphi|| <= 1 to be in the EOA.

    If phiq is in the EOA, Rphiq is obtained by linear interpolation:
        Rphiq= Rphi + A.dphi.

  3)If phiq is not in the EOA, then the mapping is computed. But as the EOA
    is a conservative approximation of the region of accuracy surrounding the
    point phi, we could expand it by comparing the computed results with the
    one obtained by linear interpolation. The error epsGrow is calculated:
        epsGrow = ||B.(dR - dRl)||,
    with dR = Rphiq - Rphi, dRl = A.dphi and B the diagonal scale factor
    matrix.
    If epsGrow <= tolerance, the EOA is too conservative and a GROW is perforned
    otherwise, the newly computed mapping is associated to the initial
    composition and added to the tree.

  4)To GROW the EOA, we expand it to include the previous EOA and the query
    point phiq. The rank-one matrix method is used. The EOA is transformed
    to a hypersphere centered at the origin. Then it is expanded to include
    the transformed point phiq' on its boundary. Then the inverse transformation
    give the modified matrix L' (see below: "Grow the EOA").


  Computation of L :
    In [1], the EOA of the constant approximation is given by
        E = {x| ||B.A/tolerance.(x-phi)|| <= 1},
    with B a scale factor diagonal matrix, A the mapping gradient matrix and
    tolerance the absolute tolerance. If we take the QR decomposition of
    (B.A)/tolerance= Q.R, with Q an orthogonal matrix and R an upper triangular
    matrix such that the EOA is described by
    (phiq-phi0)^T.R^T.R.(phiq-phi0) <= 1
    L^T = R, both Cholesky decomposition of A^T.B^T.B.A/tolerance^2
    This representation of the ellipsoid is used in [2] and in order to avoid
    large value of semi-axe length in certain direction, a Singular Value
    Decomposition (SVD) is performed on the L matrix:
        L = UDV^T,
    with the orthogonal matrix U giving the directions of the principal axes
    and 1/di the inverse of the element of the diagonal matrix D giving the
    length of the principal semi-axes. To avoid very large value of those
    length,
    di' = max(di, 1/(alphaEOA*sqrt(tolerance))), with alphaEOA = 0.1 (see [2])
    di' = max(di, 1/2), see [1]. The latter will be used in this implementation.
    And L' = UD'V^T, with D' the diagonal matrix with the modified di'.

  Grow the EOA :
    More details about the minimum-volume ellipsoid covering an ellispoid E and
    a point p are found in [3]. Here is the main steps to obtain the modified
    matrix L' describind the new ellipsoid.
        1) calculate the point p' in the transformed space :
            p' = L^T.(p-phi)
        2) compute the rank-one decomposition:
            G = I + gamma.p'.p'^T,
           with gamma = (1/|p'|-1)*1/|p'|^2
        3) compute L':
            L' = L.G.

    References:
    \verbatim
        [1] Pope, S. B. (1997).
        Computationally efficient implementation of combustion chemistry using
        in situ adaptive tabulation.
        Combustion Theory and Modelling, 1, 41-63.

        [2] Lu, L., & Pope, S. B. (2009).
        An improved algorithm for in situ adaptive tabulation.
        Journal of Computational Physics, 228(2), 361-386.

        [3] Pope, S. B. (2008).
        Algorithms for ellipsoids.
        Cornell University Report No. FDA, 08-01.
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef chemPointISAT_H
#define chemPointISAT_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
class binaryNode;

template<class CompType, class ThermoType>
class TDACChemistryModel;


/*---------------------------------------------------------------------------*\
                       Class chemPointISAT Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class chemPointISAT
{
    // Private data

        //- Pointer to the chemistryModel object
        TDACChemistryModel<CompType, ThermoType>& chemistry_;

        //- Vector storing the composition, temperature and pressure
        //  and deltaT if a variable time step is set on
        scalarField phi_;

        //- Vector storing the mapping of the composition phi
        scalarField Rphi_;

        //- LT the transpose of the L matrix describing the Ellipsoid Of
        //  Accuracy use List of Lists to be able to change size if DAC is used
        scalarSquareMatrix LT_;

        //- A the mapping gradient matrix
        scalarSquareMatrix A_;

        //- Vector storing the scale factor
        scalarField scaleFactor_;

        //- Reference to the node in the binary tree holding this chemPoint
        binaryNode<CompType, ThermoType>* node_;

        //- The size of the composition space (size of the vector phi)
        label completeSpaceSize_;

        //- Number of times the element has been grown
        label nGrowth_;

        //- Tolerance for the Ellipsoid of accuracy
        static scalar tolerance_;

        //- Number of active species stored in the chemPoint
        label nActiveSpecies_;

        //- Vectors that store the index conversion between the simplified
        //  and the complete chemical mechanism
        List<label> simplifiedToCompleteIndex_;

        label timeTag_;
        label lastTimeUsed_;

        bool toRemove_;

        label maxNumNewDim_;

        Switch printProportion_;

        //- Variable to store the number of retrieves the chemPoint
        //  will generate at each time step
        label numRetrieve_;

        //- Variable to store the number of time steps the chempoint is allowed
        //   to still live according to the maxChPLifeTime_ parameter
        label nLifeTime_;

        List<label> completeToSimplifiedIndex_;

        //- Number of equations in addition to the species eqs.
        label nAdditionalEqns_;

        label idT_;
        label idp_;
        label iddeltaT_;

        //- QR decomposition of a matrix
        //  Input : nCols cols number
        //  R the matrix to decompose
        //  QT an empty matrix that stores the transpose of the Q matrix
        //  R is returned in the given R matrix
        //  which is used to store the ellipsoid of accuracy
        void qrDecompose
        (
            const label nCols,
            scalarSquareMatrix& R
        );

        //- QR update of the matrix A
        void qrUpdate
        (
            scalarSquareMatrix& R,
            const label n,
            const scalarField& u,
            const scalarField& v
        );

        void rotate
        (
            scalarSquareMatrix& R,
            const label i,
            const scalar a,
            const scalar b,
            label n
        );


public:

    // Constructors

        //- Construct from components
        chemPointISAT
        (
            TDACChemistryModel<CompType, ThermoType>& chemistry,
            const scalarField& phi,
            const scalarField& Rphi,
            const scalarSquareMatrix& A,
            const scalarField& scaleFactor,
            const scalar& tolerance,
            const label& completeSpaceSize,
            const dictionary& coeffsDict,
            binaryNode<CompType, ThermoType>* node = nullptr
        );

        //- Construct from another chemPoint and reference to a binary node
        chemPointISAT
        (
            const chemPointISAT<CompType, ThermoType>& p,
            binaryNode<CompType, ThermoType>* node
        );

        //- Construct from another chemPoint
        chemPointISAT
        (
            chemPointISAT<CompType, ThermoType>& p
        );


    // Member functions

        //- Access to the TDACChemistryModel
        inline TDACChemistryModel<CompType, ThermoType>& chemistry()
        {
            return chemistry_;
        }

        inline label nGrowth()
        {
            return nGrowth_;
        }

        inline label& completeSpaceSize()
        {
            return completeSpaceSize_;
        }

        inline const scalarField& phi() const
        {
            return phi_;
        }

        inline const scalarField& Rphi() const
        {
            return Rphi_;
        }

        inline const scalarField& scaleFactor()
        {
            return scaleFactor_;
        }

        inline const scalar& tolerance()
        {
            return tolerance_;
        }

        static void changeTolerance(scalar newTol)
        {
            tolerance_ = newTol;
        }

        inline binaryNode<CompType, ThermoType>*& node()
        {
            return node_;
        }

        inline const scalarSquareMatrix& A() const
        {
            return A_;
        }

        inline scalarSquareMatrix& A()
        {
            return A_;
        }

        inline const scalarSquareMatrix& LT() const
        {
            return LT_;
        }

        inline scalarSquareMatrix& LT()
        {
            return LT_;
        }

        inline label nActiveSpecies()
        {
            return nActiveSpecies_;
        }

        inline List<label>& completeToSimplifiedIndex()
        {
            return completeToSimplifiedIndex_;
        }

        inline List<label>& simplifiedToCompleteIndex()
        {
            return simplifiedToCompleteIndex_;
        }

        //- Increases the number of retrieves the chempoint has generated
        void increaseNumRetrieve();

        //- Resets the number of retrieves at each time step
        void resetNumRetrieve();

        //- Increases the "counter" of the chP life
        void increaseNLifeTime();

        label simplifiedToCompleteIndex(const label i);

        inline const label& timeTag()
        {
            return timeTag_;
        }

        inline label& lastTimeUsed()
        {
            return lastTimeUsed_;
        }

        inline bool& toRemove()
        {
            return toRemove_;
        }

        inline label& maxNumNewDim()
        {
            return maxNumNewDim_;
        }

        inline const label& numRetrieve()
        {
            return numRetrieve_;
        }

        inline const label& nLifeTime()
        {
            return nLifeTime_;
        }

        inline bool variableTimeStep() const
        {
            return chemistry_.variableTimeStep();
        }

        // ISAT functions

            //- To RETRIEVE the mapping from the stored chemPoint phi, the query
            // point phiq has to be in the EOA of phi.
            // To test if phiq is in the ellipsoid:
            // ||L^T.dphi|| <= 1
            bool inEOA(const scalarField& phiq);

            //- More details about the minimum-volume ellipsoid covering an
            //  ellispoid E and a point p are found in [1].
            //  Here is the main steps to obtain the
            //  modified matrix L' describind the new ellipsoid.
            //  1) calculate the point p' in the transformed space :
            //  p' = L^T.(p-phi)
            //  2) compute the rank-one decomposition:
            //  G = I + gamma.p'.p'^T,
            //  with gamma = (1/|p'|-1)*1/|p'|^2
            //  3) compute L':
            //  L'L'^T = (L.G)(L.G)^T,
            //  L'^T is then obtained by QR decomposition of (L.G)^T = G^T.L^T
            //  [1] Stephen B. Pope, "Algorithms for ellipsoids", FDA 08-01,
            //  Cornell University, 2008
            bool grow(const scalarField& phiq);

            //- If phiq is not in the EOA, then the mapping is computed.
            //  But as the EOA is a conservative approximation of the region of
            //  accuracy surrounding the point phi, we could expand it by
            //  comparing thecomputed results with the one obtained by linear
            //  interpolation.  The error eps is calculated:
            //  eps = ||B.(dR - dRl)||,
            //  with dR = Rphiq - Rphi, dRl = A.dphi and B the diagonal scale
            //  factor matrix.
            //  If eps <= tolerance, the EOA is too conservative and a GROW is
            //  performed,
            //  otherwise, the newly computed mapping is associated to the
            //  initial composition and added to the tree.
            bool checkSolution
            (
                const scalarField& phiq,
                const scalarField& Rphiq
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "chemPointISAT.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
