/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::binaryNode<CompType, ThermoType>::binaryNode
(
)
:
    leafLeft_(nullptr),
    leafRight_(nullptr),
    nodeLeft_(nullptr),
    nodeRight_(nullptr),
    parent_(nullptr)
{}


template<class CompType, class ThermoType>
Foam::binaryNode<CompType, ThermoType>::binaryNode
(
    chemPointISAT<CompType, ThermoType>* elementLeft,
    chemPointISAT<CompType, ThermoType>* elementRight,
    binaryNode<CompType, ThermoType>* parent
)
:
    leafLeft_(elementLeft),
    leafRight_(elementRight),
    nodeLeft_(nullptr),
    nodeRight_(nullptr),
    parent_(parent),
    v_(elementLeft->completeSpaceSize(),0.0)
{
    calcV(elementLeft, elementRight, v_);
    a_ = calcA(elementLeft, elementRight);
}

template<class CompType, class ThermoType>
Foam::binaryNode<CompType, ThermoType>::binaryNode
(
    binaryNode<CompType, ThermoType> *bn
)
:
    leafLeft_(bn->leafLeft()),
    leafRight_(bn->leafRight()),
    nodeLeft_(bn->nodeLeft()),
    nodeRight_(bn->nodeRight()),
    parent_(bn->parent()),
    v_(bn->v()),
    a_(bn->a())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class CompType, class ThermoType>
void
Foam::binaryNode<CompType, ThermoType>::calcV
(
    chemPointISAT<CompType, ThermoType>*& elementLeft,
    chemPointISAT<CompType, ThermoType>*& elementRight,
    scalarField& v
)
{
    // LT is the transpose of the L matrix
    scalarSquareMatrix& LT = elementLeft->LT();
    bool mechReductionActive = elementLeft->chemistry().mechRed()->active();
    // difference of composition in the full species domain
    scalarField phiDif(elementRight->phi() - elementLeft->phi());
    const scalarField& scaleFactor(elementLeft->scaleFactor());
    scalar epsTol = elementLeft->tolerance();

    // v = LT.T()*LT*phiDif
    for (label i=0; i<elementLeft->completeSpaceSize(); i++)
    {
        label si = i;
        bool outOfIndexI = true;
        if (mechReductionActive)
        {
            if (i<elementLeft->completeSpaceSize()-2)
            {
                si = elementLeft->completeToSimplifiedIndex()[i];
                outOfIndexI = (si==-1);
            }
            else// temperature and pressure
            {
                outOfIndexI = false;
                label dif = i-(elementLeft->completeSpaceSize()-2);
                si = elementLeft->nActiveSpecies()+dif;
            }
        }
        if (!mechReductionActive || (mechReductionActive && !(outOfIndexI)))
        {
            v[i]=0.0;
            for (label j=0; j<elementLeft->completeSpaceSize(); j++)
            {
                label sj = j;
                bool outOfIndexJ = true;
                if (mechReductionActive)
                {
                    if (j<elementLeft->completeSpaceSize()-2)
                    {
                        sj = elementLeft->completeToSimplifiedIndex()[j];
                        outOfIndexJ = (sj==-1);
                    }
                    else
                    {
                        outOfIndexJ = false;
                        label dif = j-(elementLeft->completeSpaceSize()-2);
                        sj = elementLeft->nActiveSpecies()+dif;
                    }
                }
                if
                (
                    !mechReductionActive
                  ||(mechReductionActive && !(outOfIndexJ))
                )
                {
                    // since L is a lower triangular matrix k=0->min(i, j)
                    for (label k=0; k<=min(si, sj); k++)
                    {
                        v[i] += LT(k, si)*LT(k, sj)*phiDif[j];
                    }
                }
            }
        }
        else
        {
            // when it is an inactive species the diagonal element of LT is
            //  1/(scaleFactor*epsTol)
            v[i] = phiDif[i]/sqr(scaleFactor[i]*epsTol);
        }
    }
}


template<class CompType, class ThermoType>
Foam::scalar Foam::binaryNode<CompType, ThermoType>::calcA
(
    chemPointISAT<CompType, ThermoType>* elementLeft,
    chemPointISAT<CompType, ThermoType>* elementRight
)
{
    scalar a = 0.0;
    scalarField phih((elementLeft->phi()+elementRight->phi())/2);
    label completeSpaceSize = elementLeft->completeSpaceSize();
    for (label i=0; i<completeSpaceSize; i++)
    {
        a += v_[i]*phih[i];
    }
    return a;
}


// ************************************************************************* //
