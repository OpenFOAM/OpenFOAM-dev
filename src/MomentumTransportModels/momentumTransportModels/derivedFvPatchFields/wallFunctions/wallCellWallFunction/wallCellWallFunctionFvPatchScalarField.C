/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "wallCellWallFunctionFvPatchScalarField.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::wallCellWallFunctionFvPatchScalarField::tol_ = 1e-1;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::wallCellWallFunctionFvPatchScalarField::initMaster()
{
    volScalarField& vf =
        const_cast<volScalarField&>
        (
            static_cast<const volScalarField&>
            (
                internalField()
            )
        );

    volScalarField::Boundary& bf = vf.boundaryFieldRef();

    // Choose a master patch if one has not previously been chosen
    if (masterPatchi_ == -1)
    {
        label masterPatchi = -1;

        forAll(bf, patchi)
        {
            if (isA<wallCellWallFunctionFvPatchScalarField>(bf[patchi]))
            {
                wallCellWallFunctionFvPatchScalarField& epsilonPf =
                    refCast<wallCellWallFunctionFvPatchScalarField>(bf[patchi]);

                if (masterPatchi == -1)
                {
                    masterPatchi = patchi;
                }

                epsilonPf.masterPatchi_ = masterPatchi;
            }
        }
    }

    // Additional initialisation is only necessary on the master patch
    if (patch().index() != masterPatchi_)
    {
        return;
    }

    // Initialise at the start, and re-initialise following any mesh changes
    if (wallCellsPtr_.valid() && !vf.mesh().changing())
    {
        return;
    }

    // Construct flag lists for whether a cell or a (poly) boundary face is
    // associated with a wall function
    boolList cellIsWall(vf.mesh().nCells(), false);
    boolList bFaceIsWall
    (
        vf.mesh().nFaces() - vf.mesh().nInternalFaces(),
        false
    );
    forAll(bf, patchi)
    {
        if (isA<wallCellWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            UIndirectList<bool>
            (
                cellIsWall,
                vf.mesh().boundary()[patchi].faceCells()
            ) = true;
            UIndirectList<bool>
            (
                bFaceIsWall,
                vf.mesh().polyFacesBf()[patchi] - vf.mesh().nInternalFaces()
            ) = true;
        }
    }

    // Identify the wall cells
    wallCellsPtr_.reset(new labelList(findIndices(cellIsWall, true)));

    // Construct the poly wall areas for each cell
    scalarField cellMagWallArea(vf.mesh().nCells(), 0);
    forAll(bFaceIsWall, bFacei)
    {
        if (bFaceIsWall[bFacei])
        {
            const label facei = bFacei + vf.mesh().nInternalFaces();
            const label celli = vf.mesh().faceOwner()[facei];
            cellMagWallArea[celli] += vf.mesh().magFaceAreas()[facei];
        }
    }

    // Construct the fv wall areas for each cell
    scalarField cellMagSw(vf.mesh().nCells(), 0);
    forAll(bf, patchi)
    {
        const fvPatch& fvp = vf.mesh().boundary()[patchi];

        if (isA<wallCellWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            patchFieldAddToCellField(fvp, fvp.magSf(), cellMagSw);
        }
    }

    // Divide the areas to get the fraction of the cell values that is to be
    // constrained. Modify this by the tolerance so that small fractions do not
    // constrain anything.
    tmp<scalarField> tWallCellFraction
    (
        scalarField(cellMagSw, wallCellsPtr_())
       /scalarField(cellMagWallArea, wallCellsPtr_())
    );
    wallCellFractionPtr_.reset
    (
        max((tWallCellFraction - tol_)/(1 - tol_), scalar(0)).ptr()
    );
}


void Foam::wallCellWallFunctionFvPatchScalarField::patchFieldAddToCellField
(
    const fvPatch& fvp,
    const scalarField& pf,
    scalarField& vf
)
{
    const labelUList& patchFaceCells = fvp.faceCells();

    forAll(patchFaceCells, patchFacei)
    {
        vf[patchFaceCells[patchFacei]] += pf[patchFacei];
    }
}


Foam::tmp<Foam::scalarField>
Foam::wallCellWallFunctionFvPatchScalarField::patchFieldsToWallCellField
(
    const PtrList<scalarField>& pfs
) const
{
    const volScalarField& vf =
        static_cast<const volScalarField&>(internalField());

    const volScalarField::Boundary& bf = vf.boundaryField();

    scalarField cellMagSw(vf.mesh().nCells(), scalar(0));
    scalarField cellMagSwValue(vf.mesh().nCells(), scalar(0));

    forAll(bf, patchi)
    {
        const fvPatch& fvp = vf.mesh().boundary()[patchi];

        if (isA<wallCellWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            patchFieldAddToCellField
            (
                fvp,
                fvp.magSf(),
                cellMagSw
            );
            patchFieldAddToCellField
            (
                fvp,
                fvp.magSf()*pfs[patchi],
                cellMagSwValue
            );
        }
    }

    return
        scalarField(cellMagSwValue, wallCellsPtr_())
       /max(scalarField(cellMagSw, wallCellsPtr_()), vSmall);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallCellWallFunctionFvPatchScalarField::
wallCellWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict, valueRequired),
    masterPatchi_(-1),
    wallCellsPtr_(nullptr),
    wallCellFractionPtr_(nullptr)
{}


Foam::wallCellWallFunctionFvPatchScalarField::
wallCellWallFunctionFvPatchScalarField
(
    const wallCellWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper,
    const bool mappingRequired
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper, mappingRequired),
    masterPatchi_(-1),
    wallCellsPtr_(nullptr),
    wallCellFractionPtr_(nullptr)
{}


Foam::wallCellWallFunctionFvPatchScalarField::
wallCellWallFunctionFvPatchScalarField
(
    const wallCellWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ewfpsf, iF),
    masterPatchi_(-1),
    wallCellsPtr_(nullptr),
    wallCellFractionPtr_(nullptr)
{}


// ************************************************************************* //
