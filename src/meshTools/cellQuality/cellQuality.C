/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "cellQuality.H"
#include "unitConversion.H"
#include "SubField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellQuality::cellQuality(const polyMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::cellQuality::nonOrthogonality() const
{
    tmp<scalarField> tresult
    (
        new scalarField
        (
            mesh_.nCells(), 0.0
        )
    );

    scalarField& result = tresult();

    scalarField sumArea(mesh_.nCells(), 0.0);

    const vectorField& centres = mesh_.cellCentres();
    const vectorField& areas = mesh_.faceAreas();

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    forAll(nei, faceI)
    {
        vector d = centres[nei[faceI]] - centres[own[faceI]];
        vector s = areas[faceI];
        scalar magS = mag(s);

        scalar cosDDotS =
            radToDeg(Foam::acos(min(1.0, (d & s)/(mag(d)*magS + VSMALL))));

        result[own[faceI]] = max(cosDDotS, result[own[faceI]]);

        result[nei[faceI]] = max(cosDDotS, result[nei[faceI]]);
    }

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const labelUList& faceCells =
            mesh_.boundaryMesh()[patchI].faceCells();

        const vectorField::subField faceCentres =
            mesh_.boundaryMesh()[patchI].faceCentres();

        const vectorField::subField faceAreas =
            mesh_.boundaryMesh()[patchI].faceAreas();

        forAll(faceCentres, faceI)
        {
            vector d = faceCentres[faceI] - centres[faceCells[faceI]];
            vector s = faceAreas[faceI];
            scalar magS = mag(s);

            scalar cosDDotS =
                radToDeg(Foam::acos(min(1.0, (d & s)/(mag(d)*magS + VSMALL))));

            result[faceCells[faceI]] = max(cosDDotS, result[faceCells[faceI]]);
        }
    }

    return tresult;
}


Foam::tmp<Foam::scalarField> Foam::cellQuality::skewness() const
{
    tmp<scalarField> tresult
    (
        new scalarField
        (
            mesh_.nCells(), 0.0
        )
    );
    scalarField& result = tresult();

    scalarField sumArea(mesh_.nCells(), 0.0);

    const vectorField& cellCtrs = mesh_.cellCentres();
    const vectorField& faceCtrs = mesh_.faceCentres();
    const vectorField& areas = mesh_.faceAreas();

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    forAll(nei, faceI)
    {
        scalar dOwn = mag
        (
            (faceCtrs[faceI] - cellCtrs[own[faceI]]) & areas[faceI]
        )/mag(areas[faceI]);

        scalar dNei = mag
        (
            (cellCtrs[nei[faceI]] - faceCtrs[faceI]) & areas[faceI]
        )/mag(areas[faceI]);

        point faceIntersection =
            cellCtrs[own[faceI]]
          + (dOwn/(dOwn+dNei))*(cellCtrs[nei[faceI]] - cellCtrs[own[faceI]]);

        scalar skewness =
            mag(faceCtrs[faceI] - faceIntersection)
           /(mag(cellCtrs[nei[faceI]] - cellCtrs[own[faceI]]) + VSMALL);

        result[own[faceI]] = max(skewness, result[own[faceI]]);

        result[nei[faceI]] = max(skewness, result[nei[faceI]]);
    }

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const labelUList& faceCells =
            mesh_.boundaryMesh()[patchI].faceCells();

        const vectorField::subField faceCentres =
            mesh_.boundaryMesh()[patchI].faceCentres();

        const vectorField::subField faceAreas =
            mesh_.boundaryMesh()[patchI].faceAreas();

        forAll(faceCentres, faceI)
        {
            vector n = faceAreas[faceI]/mag(faceAreas[faceI]);

            point faceIntersection =
                cellCtrs[faceCells[faceI]]
              + ((faceCentres[faceI] - cellCtrs[faceCells[faceI]])&n)*n;

            scalar skewness =
                mag(faceCentres[faceI] - faceIntersection)
               /(
                    mag(faceCentres[faceI] - cellCtrs[faceCells[faceI]])
                  + VSMALL
                );

            result[faceCells[faceI]] = max(skewness, result[faceCells[faceI]]);
        }
    }

    return tresult;
}


Foam::tmp<Foam::scalarField> Foam::cellQuality::faceNonOrthogonality() const
{
    tmp<scalarField> tresult
    (
        new scalarField
        (
            mesh_.nFaces(), 0.0
        )
    );
    scalarField& result = tresult();


    const vectorField& centres = mesh_.cellCentres();
    const vectorField& areas = mesh_.faceAreas();

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    forAll(nei, faceI)
    {
        vector d = centres[nei[faceI]] - centres[own[faceI]];
        vector s = areas[faceI];
        scalar magS = mag(s);

        scalar cosDDotS =
            radToDeg(Foam::acos(min(1.0, (d & s)/(mag(d)*magS + VSMALL))));

        result[faceI] = cosDDotS;
    }

    label globalFaceI = mesh_.nInternalFaces();

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const labelUList& faceCells =
            mesh_.boundaryMesh()[patchI].faceCells();

        const vectorField::subField faceCentres =
            mesh_.boundaryMesh()[patchI].faceCentres();

        const vectorField::subField faceAreas =
            mesh_.boundaryMesh()[patchI].faceAreas();

        forAll(faceCentres, faceI)
        {
            vector d = faceCentres[faceI] - centres[faceCells[faceI]];
            vector s = faceAreas[faceI];
            scalar magS = mag(s);

            scalar cosDDotS =
                radToDeg(Foam::acos(min(1.0, (d & s)/(mag(d)*magS + VSMALL))));

            result[globalFaceI++] = cosDDotS;
        }
    }

    return tresult;
}


Foam::tmp<Foam::scalarField> Foam::cellQuality::faceSkewness() const
{
    tmp<scalarField> tresult
    (
        new scalarField
        (
            mesh_.nFaces(), 0.0
        )
    );
    scalarField& result = tresult();


    const vectorField& cellCtrs = mesh_.cellCentres();
    const vectorField& faceCtrs = mesh_.faceCentres();
    const vectorField& areas = mesh_.faceAreas();

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    forAll(nei, faceI)
    {
        scalar dOwn = mag
        (
            (faceCtrs[faceI] - cellCtrs[own[faceI]]) & areas[faceI]
        )/mag(areas[faceI]);

        scalar dNei = mag
        (
            (cellCtrs[nei[faceI]] - faceCtrs[faceI]) & areas[faceI]
        )/mag(areas[faceI]);

        point faceIntersection =
            cellCtrs[own[faceI]]
          + (dOwn/(dOwn+dNei))*(cellCtrs[nei[faceI]] - cellCtrs[own[faceI]]);

        result[faceI] =
            mag(faceCtrs[faceI] - faceIntersection)
           /(mag(cellCtrs[nei[faceI]] - cellCtrs[own[faceI]]) + VSMALL);
    }


    label globalFaceI = mesh_.nInternalFaces();

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const labelUList& faceCells =
            mesh_.boundaryMesh()[patchI].faceCells();

        const vectorField::subField faceCentres =
            mesh_.boundaryMesh()[patchI].faceCentres();

        const vectorField::subField faceAreas =
            mesh_.boundaryMesh()[patchI].faceAreas();

        forAll(faceCentres, faceI)
        {
            vector n = faceAreas[faceI]/mag(faceAreas[faceI]);

            point faceIntersection =
                cellCtrs[faceCells[faceI]]
              + ((faceCentres[faceI] - cellCtrs[faceCells[faceI]])&n)*n;

            result[globalFaceI++] =
                mag(faceCentres[faceI] - faceIntersection)
               /(
                    mag(faceCentres[faceI] - cellCtrs[faceCells[faceI]])
                  + VSMALL
                );
        }
    }

    return tresult;
}


// ************************************************************************* //
