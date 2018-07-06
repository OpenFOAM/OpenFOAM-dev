/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "fileControl.H"
#include "addToRunTimeSelectionTable.H"
#include "tetPointRef.H"
#include "scalarList.H"
#include "vectorTools.H"
#include "pointIOField.H"
#include "scalarIOField.H"
#include "triadIOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(fileControl, 0);
addToRunTimeSelectionTable
(
    cellSizeAndAlignmentControl,
    fileControl,
    dictionary
);

}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileControl::fileControl
(
    const Time& runTime,
    const word& name,
    const dictionary& controlFunctionDict,
    const conformationSurfaces& geometryToConformTo,
    const scalar& defaultCellSize
)
:
    cellSizeAndAlignmentControl
    (
        runTime,
        name,
        controlFunctionDict,
        geometryToConformTo,
        defaultCellSize
    ),
    pointsFile_(controlFunctionDict.lookup("pointsFile")),
    sizesFile_(controlFunctionDict.lookup("sizesFile")),
    alignmentsFile_(controlFunctionDict.lookup("alignmentsFile")),
    maxPriority_(readLabel(controlFunctionDict.lookup("priority")))
{
    Info<< indent << "Loading " << name << " from file:" << nl
        << indent << "    priority   : " << maxPriority_ << nl
        << indent << "    points     : " << pointsFile_ << nl
        << indent << "    sizes      : " << sizesFile_ << nl
        << indent << "    alignments : " << alignmentsFile_
        << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileControl::~fileControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
//
//Foam::scalar Foam::fileControl::cellSize(const point& pt) const
//{
//    FixedList<scalar, 4> bary;
//    Cell_handle ch;
//
//    triangulatedMesh_.barycentricCoords(pt, bary, ch);
//
//    scalar size = 0;
//    forAll(bary, pI)
//    {
//        size += bary[pI]*ch->vertex(pI)->size();
//    }
//
//    return size;
//}
//
//
////- Return the cell alignment at the given location
//Foam::tensor Foam::fileControl::cellAlignment(const point& pt) const
//{
//    FixedList<scalar, 4> bary;
//    Cell_handle ch;
//
//    triangulatedMesh_.barycentricCoords(pt, bary, ch);
//
//    label nearest = 0;
//
//    tensor alignment = Zero;
//    forAll(bary, pI)
//    {
//        // alignment += bary[pI]*ch->vertex(pI)->alignment();
//
//        // Find nearest point
//        if (bary[pI] > nearest)
//        {
//            alignment = ch->vertex(pI)->alignment();
//            nearest = bary[pI];
//        }
//    }
//
//    return alignment;
//}
//
//
////- Return the cell alignment at the given location
//void Foam::fileControl::cellSizeAndAlignment
//(
//    const point& pt,
//    scalar& size,
//    tensor& alignment
//) const
//{
//    FixedList<scalar, 4> bary;
//    Cell_handle ch;
//
//    triangulatedMesh_.barycentricCoords(pt, bary, ch);
//
//    size = 0;
//    forAll(bary, pI)
//    {
//        size += bary[pI]*ch->vertex(pI)->size();
//    }
//
////    alignment = Zero;
////    forAll(bary, pI)
////    {
////        alignment += bary[pI]*ch->vertex(pI)->alignment();
////    }
//
//    alignment = cellAlignment(pt);
//}


void Foam::fileControl::cellSizeFunctionVertices
(
    DynamicList<Foam::point>& pts,
    DynamicList<scalar>& sizes
) const
{
    return;
}


void Foam::fileControl::initialVertices
(
    pointField& pts,
    scalarField& sizes,
    triadField& alignments
) const
{
    Info<< "    Reading points from file     : " << pointsFile_ << endl;

    pointIOField pointsTmp
    (
        IOobject
        (
            pointsFile_,
            runTime_.constant(),
            runTime_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    pts.transfer(pointsTmp);

    Info<< "    Reading sizes from file      : " << sizesFile_ << endl;

    scalarIOField sizesTmp
    (
        IOobject
        (
            sizesFile_,
            runTime_.constant(),
            runTime_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    sizes.transfer(sizesTmp);

    Info<< "    Reading alignments from file : " << alignmentsFile_ << endl;

    triadIOField alignmentsTmp
    (
        IOobject
        (
            alignmentsFile_,
            runTime_.constant(),
            runTime_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    alignments.transfer(alignmentsTmp);

    if ((pts.size() != sizes.size()) || (pts.size() != alignments.size()))
    {
        FatalErrorInFunction
            << "Size of list of points, sizes and alignments do not match:"
            << nl
            << "Points size     = " << pts.size() << nl
            << "Sizes size      = " << sizes.size() << nl
            << "Alignments size = " << alignments.size()
            << abort(FatalError);
    }
}


// ************************************************************************* //
