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

#include "binaryNode.H"
#include "ISAT.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binaryNode::binaryNode()
:
    leafLeft_(nullptr),
    leafRight_(nullptr),
    nodeLeft_(nullptr),
    nodeRight_(nullptr),
    parent_(nullptr)
{}


Foam::binaryNode::binaryNode
(
    chemPointISAT* elementLeft,
    chemPointISAT* elementRight,
    binaryNode* parent
)
:
    leafLeft_(elementLeft),
    leafRight_(elementRight),
    nodeLeft_(nullptr),
    nodeRight_(nullptr),
    parent_(parent),
    v_(elementLeft->completeSpaceSize(), 0)
{
    calcV(*elementLeft, *elementRight, v_);
    a_ = calcA(*elementLeft, *elementRight);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::binaryNode::calcV
(
    const chemPointISAT& elementLeft,
    const chemPointISAT& elementRight,
    scalarField& v
)
{
    // LT is the transpose of the L matrix
    const scalarSquareMatrix& LT = elementLeft.LT();
    const bool reduction = elementLeft.table().reduction();

    // Difference of composition in the full species domain
    const scalarField& phil(elementLeft.phi());
    const scalarField& phir(elementRight.phi());
    const scalarField& scaleFactor(elementLeft.scaleFactor());
    const scalar epsTol = elementLeft.tolerance();

    // v = LT.T()*LT*(phir - phil)
    for (label i=0; i<elementLeft.completeSpaceSize(); i++)
    {
        label si = i;
        bool outOfIndexI = true;
        if (reduction)
        {
            if (i<elementLeft.completeSpaceSize() - 3)
            {
                si = elementLeft.completeToSimplifiedIndex()[i];
                outOfIndexI = (si == -1);
            }
            else // temperature and pressure
            {
                outOfIndexI = false;
                const label dif = i - (elementLeft.completeSpaceSize() - 3);
                si = elementLeft.nActive() + dif;
            }
        }
        if (!reduction || (reduction && !(outOfIndexI)))
        {
            v[i] = 0;
            for (label j=0; j<elementLeft.completeSpaceSize(); j++)
            {
                label sj = j;
                bool outOfIndexJ = true;

                if (reduction)
                {
                    if (j < elementLeft.completeSpaceSize() - 3)
                    {
                        sj = elementLeft.completeToSimplifiedIndex()[j];
                        outOfIndexJ = (sj==-1);
                    }
                    else
                    {
                        outOfIndexJ = false;
                        const label dif =
                            j - (elementLeft.completeSpaceSize() - 3);
                        sj = elementLeft.nActive() + dif;
                    }
                }

                if (!reduction ||(reduction && !(outOfIndexJ)))
                {
                    // Since L is a lower triangular matrix k=0->min(i, j)
                    for (label k=0; k<=min(si, sj); k++)
                    {
                        v[i] += LT(k, si)*LT(k, sj)*(phir[j] - phil[j]);
                    }
                }
            }
        }
        else
        {
            // When it is an inactive species the diagonal element of LT is
            // 1/(scaleFactor*epsTol)
            v[i] = (phir[i] - phil[i])/sqr(scaleFactor[i]*epsTol);
        }
    }
}


// ************************************************************************* //
