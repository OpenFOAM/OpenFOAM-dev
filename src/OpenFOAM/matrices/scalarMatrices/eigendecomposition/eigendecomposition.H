/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::eigendecomposition

Description
    Eigen decomposition of a square matrix.

    Calculates the eigenvalues and eigenvectors of a square scalar matrix
    matrix.

    If the matrix \f$ A \f$ is symmetric, then \f$ A = VDV^-1 \f$ where the
    eigenvalue matrix \f$ D \f$ is diagonal and the eigenvector matrix \f$ V \f$
    is orthogonal. That is, the diagonal values of \f$ D \f$ are the
    eigenvalues, and \f$ VV^-1 = I \f$, where \f$ I \f$ is the identity matrix.
    The columns of V represent the eigenvectors in the sense that \f$ AV = VD
    \f$.

    If \f$ A \f$ is not symmetric, then the eigenvalue matrix \f$ D \f$ is
    block diagonal with the real eigenvalues in 1-by-1 blocks and any complex
    eigenvalues, \f$ a + ib \f$, in 2-by-2 blocks, \f$ [a, b; -b, a] \f$.  That
    is, if the complex eigenvalues look like
    \verbatim
          u + iv     .        .          .      .    .
            .      u - iv     .          .      .    .
            .        .      a + ib       .      .    .
            .        .        .        a - ib   .    .
            .        .        .          .      x    .
            .        .        .          .      .    y
    \endverbatim
    then \f$ D \f$ looks like
    \verbatim
            u        v        .          .      .    .
           -v        u        .          .      .    .
            .        .        a          b      .    .
            .        .       -b          a      .    .
            .        .        .          .      x    .
            .        .        .          .      .    y
    \endverbatim
    This keeps \f$ V \f$ a real matrix in both symmetric and non-symmetric
    cases, and \f$ AV = VD \f$.

    The matrix \f$ V \f$ may be badly conditioned, or even singular, so the
    validity of the equation \f$ A = VDV^-1 \f$ depends upon the condition
    number of \f$ V \f$.

    Adapted from the TNT C++ implementation (http://math.nist.gov/tnt) of the
    JAMA Java implementation (http://math.nist.gov/javanumerics/jama) of the
    algorithms in:
    \verbatim
        Wilkinson, J. H., & Reinsch, C. (1971).
        Handbook for Automatic Computation: Volume II:
        Linear Algebra (Vol. 186). Springer-Verlag.
    \endverbatim

SourceFiles
    eigendecomposition.C

\*---------------------------------------------------------------------------*/

#ifndef eigendecomposition_H
#define eigendecomposition_H

#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                      Class eigendecomposition Declaration
\*---------------------------------------------------------------------------*/

class eigendecomposition
{
    // Private Data

        //- Is the matrix A symmetric
        bool symmetric_;

        //- Real part of the eigenvalues
        scalarField d_;

        //- Imaginary part of the eigenvalues
        scalarField e_;

        //- Matrix for eigenvectors
        scalarSquareMatrix V_;

        //- Matrix for the asymmetric Hessenberg form
        //  Allocated if required
        scalarSquareMatrix H_;

        //- Working storage for asymmetric algorithm
        //  Allocated if required
        scalarField ort_;


    // Private Member Functions

        //- Symmetric Householder reduction to tridiagonal form
        void tred2();

        //- Symmetric tridiagonal QL algorithm
        void tql2();

        //- Asymmetric reduction to Hessenberg form
        void orthes();

        //- Complex scalar division
        inline void cdiv
        (
            scalar& cdivr,
            scalar& cdivi,
            const scalar xr,
            const scalar xi,
            const scalar yr,
            const scalar yi
        );

        //- Asymmetric reduction from Hessenberg to real Schur form.
        void hqr2();


public:

    // Constructors

        //- Construct the eigen decomposition of matrix A
        eigendecomposition(const scalarSquareMatrix& A);


    // Member Functions

        //- Return the real part of the eigenvalues
        const scalarField& d() const
        {
            return d_;
        }

        //- Return the imaginary part of the eigenvalues
        const scalarField& e() const
        {
            return e_;
        }

        //- Return the block diagonal eigenvalue matrix in D
        void D(scalarSquareMatrix& D) const;

        //- Return the eigenvector matrix
        const scalarSquareMatrix& V() const
        {
            return V_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
