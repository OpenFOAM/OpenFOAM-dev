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

\*---------------------------------------------------------------------------*/

#include "QRMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MatrixType>
Foam::QRMatrix<MatrixType>::QRMatrix(const MatrixType& M)
{
    decompose(M);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class MatrixType>
void Foam::QRMatrix<MatrixType>::decompose(const MatrixType& M)
{
    const label m = M.m();
    const label n = M.n();

    // Initialize the R-matrix to M
    R_ = M;

    // Initialize the Q-matrix to I
    Q_.setSize(m);
    Q_ = I;

    // Pre-allocate temporary storage for the Householder steps
    QMatrixType Qk(m);
    QMatrixType Rk(m);
    Field<cmptType> uk(m);

    for (label k=0; k<n; k++)
    {
        // alpha = -|column k of Rk|
        cmptType alpha = Zero;
        for (label bi=k; bi<m; bi++)
        {
            alpha += sqr(R_(bi, k));
        }
        alpha = sqrt(alpha);

        if (R_(k, k) > 0)
        {
            alpha = -alpha;
        }

        // uk = column k of Rk - alpha*ek
        // rSumSqrUk = 2/sum(sqr(uk))
        uk[k] = R_(k, k) - alpha;
        cmptType rSumSqrUk = sqr(uk[k]);
        for (label bi=k+1; bi<m; bi++)
        {
            uk[bi] = R_(bi, k);
            rSumSqrUk += sqr(uk[bi]);
        }
        rSumSqrUk = 2/rSumSqrUk;

        // Qk = I - 2*u*uT/sum(sqr(uk))
        for (label bi=k; bi<m; bi++)
        {
            for (label bj=k; bj<m; bj++)
            {
                Qk(bi, bj) = -rSumSqrUk*uk[bi]*uk[bj];
            }
            Qk(bi, bi) += 1;
        }

        // Rk = Qk*R
        for (label bi=k; bi<m; bi++)
        {
            for (label bk=k; bk<n; bk++)
            {
                Rk(bi, bk) = Zero;
                for (label bj=k; bj<n; bj++)
                {
                    Rk(bi, bk) += Qk(bi, bj)*R_(bj, bk);
                }
            }
        }

        // Ensure diagonal is positive
        if (R_(k, k) < 0)
        {
            // R = -Rk
            // Q = -Q
            for (label bi=k; bi<m; bi++)
            {
                for (label bj=k; bj<n; bj++)
                {
                    R_(bi, bj) = -Rk(bi, bj);
                }
                for (label bj=k; bj<m; bj++)
                {
                    Q_(bi, bj) = -Q_(bi, bj);
                }
            }
        }
        else
        {
            // R = Rk
            for (label bi=k; bi<m; bi++)
            {
                for (label bj=k; bj<n; bj++)
                {
                    R_(bi, bj) = Rk(bi, bj);
                }
            }
        }

        // Q = Q*Qk (using Rk as temporary storage)
        for (label bi=0; bi<m; bi++)
        {
            for (label bk=k; bk<m; bk++)
            {
                Rk(bi, bk) = Zero;
                for (label bj=k; bj<m; bj++)
                {
                    Rk(bi, bk) += Q_(bi, bj)*Qk(bj, bk);
                }
            }
        }
        for (label bi=0; bi<m; bi++)
        {
            for (label bj=k; bj<n; bj++)
            {
                Q_(bi, bj) = Rk(bi, bj);
            }
        }
    }
}


template<class MatrixType>
void Foam::QRMatrix<MatrixType>::solvex
(
    Field<cmptType>& x
) const
{
    const label n = R_.n();

    for (int i=n-1; i>=0; i--)
    {
        cmptType sum = x[i];

        for (label j=i+1; j<n; j++)
        {
            sum -= x[j]*R_(i, j);
        }

        if (mag(R_(i, i)) < small)
        {
            FatalErrorInFunction
                << "Back-substitution failed due to small diagonal"
                << abort(FatalError);
        }

        x[i] = sum/R_(i, i);
    }
}


template<class MatrixType>
void Foam::QRMatrix<MatrixType>::solve
(
    Field<cmptType>& x,
    const Field<cmptType>& source
) const
{
    const label m = Q_.m();

    // x = Q_.T()*source;
    for (label i=0; i<m; i++)
    {
        x[i] = 0;
        for (label j=0; j<m; j++)
        {
            x[i] += Q_(j, i)*source[j];
        }
    }

    solvex(x);
}


template<class MatrixType>
Foam::tmp<Foam::Field<typename MatrixType::cmptType>>
Foam::QRMatrix<MatrixType>::solve
(
    const Field<cmptType>& source
) const
{
    tmp<Field<cmptType>> tx(new Field<cmptType>(Q_.m()));
    Field<cmptType>& x = tx.ref();

    solve(x, source);

    return tx;
}


template<class MatrixType>
typename Foam::QRMatrix<MatrixType>::QMatrixType
Foam::QRMatrix<MatrixType>::inv() const
{
    const label m = Q_.m();

    Field<cmptType> x(m);
    QMatrixType inv(m);

    for (label i=0; i<m; i++)
    {
        for (label j=0; j<m; j++)
        {
            x[j] = Q_(i, j);
        }
        solvex(x);
        inv.block(m, 1, 0, i) = x;
    }

    return inv;
}


// ************************************************************************* //
