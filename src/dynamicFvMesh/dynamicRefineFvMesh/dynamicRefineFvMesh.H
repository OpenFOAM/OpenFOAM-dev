/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::dynamicRefineFvMesh

Description
    Dynamic mesh refinement/unrefinement based on volScalarField values.

    Refinement can optionally be specified in a cellZone or in multiple
    regions, each controlled by a different volScalarField.

Usage
    Example of single field based refinement in all cells:
    \verbatim
        dynamicFvMesh   dynamicRefineFvMesh;

        // How often to refine
        refineInterval  1;

        // Field to be refinement on
        field           alpha.water;

        // Refine field in between lower..upper
        lowerRefineLevel 0.001;
        upperRefineLevel 0.999;

        // If value < unrefineLevel unrefine
        unrefineLevel   10;

        // If value < unrefineLevel (default=great) unrefine
        // unrefineLevel   10;
        // Have slower than 2:1 refinement
        nBufferLayers   1;

        // Refine cells only up to maxRefinement levels
        maxRefinement   2;

        // Stop refinement if maxCells reached
        maxCells        200000;

        // Flux field and corresponding velocity field. Fluxes on changed
        // faces get recalculated by interpolating the velocity. Use 'none'
        // on surfaceScalarFields that do not need to be reinterpolated.
        correctFluxes
        (
            (phi none)
            (nHatf none)
            (rhoPhi none)
            (alphaPhi.water none)
            (ghf none)
        );

        // Write the refinement level as a volScalarField
        dumpLevel       true;
    \endverbatim

    Example of single field based refinement in two regions:
    \verbatim
        dynamicFvMesh   dynamicRefineFvMesh;

        // How often to refine
        refineInterval  1;

        refinementRegions
        {
            region1
            {
                cellZone        refinementRegion1;

                // Field to be refinement on
                field           alpha.water;

                // Refine field in between lower..upper
                lowerRefineLevel 0.001;
                upperRefineLevel 0.999;

                // Refine cells only up to maxRefinement levels
                maxRefinement   1;

                // If value < unrefineLevel unrefine
                unrefineLevel   10;
            }

            region2
            {
                cellZone        refinementRegion2;

                // Field to be refinement on
                field           alpha.water;

                // Refine field in between lower..upper
                lowerRefineLevel 0.001;
                upperRefineLevel 0.999;

                // Refine cells only up to maxRefinement levels
                maxRefinement   2;

                // If value < unrefineLevel unrefine
                unrefineLevel   10;
            }
        }

        // If value < unrefineLevel (default=great) unrefine
        // unrefineLevel   10;
        // Have slower than 2:1 refinement
        nBufferLayers   1;

        // Stop refinement if maxCells reached
        maxCells        200000;

        // Flux field and corresponding velocity field. Fluxes on changed
        // faces get recalculated by interpolating the velocity. Use 'none'
        // on surfaceScalarFields that do not need to be reinterpolated.
        correctFluxes
        (
            (phi none)
            (nHatf none)
            (rhoPhi none)
            (alphaPhi.water none)
            (ghf none)
        );

        // Write the refinement level as a volScalarField
        dumpLevel       true;
    \endverbatim


SourceFiles
    dynamicRefineFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicRefineFvMesh_H
#define dynamicRefineFvMesh_H

#include "dynamicFvMesh.H"
#include "hexRef8.H"
#include "PackedBoolList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class dynamicRefineFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicRefineFvMesh
:
    public dynamicFvMesh
{
    // Private Member Data

        //- Mesh cutting engine
        hexRef8 meshCutter_;

        //- Dump cellLevel for postprocessing
        Switch dumpLevel_;

        //- Fluxes to map
        HashTable<word> correctFluxes_;

        //- Number of refinement/unrefinement steps done so far.
        label nRefinementIterations_;

        //- Protected cells (usually since not hexes)
        PackedBoolList protectedCells_;


    // Private Member Functions

        //- Count set/unset elements in packedlist.
        static label count(const PackedBoolList&, const unsigned int);

        //- Calculate cells that cannot be refined since would trigger
        //  refinement of protectedCells_ (since 2:1 refinement cascade)
        void calculateProtectedCells(PackedBoolList& unrefineableCells) const;

        //- Read the projection parameters from dictionary
        void readDict();


        //- Refine cells. Update mesh and fields.
        autoPtr<mapPolyMesh> refine(const labelList&);

        //- Unrefine cells. Gets passed in centre points of cells to combine.
        autoPtr<mapPolyMesh> unrefine(const labelList&);


        // Selection of cells to un/refine

            const cellZone& findCellZone
            (
                const word& cellZoneName
            ) const;

            scalarField cellToPoint(const scalarField& vFld) const;

            scalarField error
            (
                const scalarField& fld,
                const scalar minLevel,
                const scalar maxLevel
            ) const;

            scalarField error
            (
                const scalarField& fld,
                const labelList& cells,
                const scalar minLevel,
                const scalar maxLevel
            ) const;

            //- Select candidate cells for refinement
            virtual void selectRefineCandidates
            (
                PackedBoolList& candidateCell,
                const scalar lowerRefineLevel,
                const scalar upperRefineLevel,
                const scalar maxRefinement,
                const scalarField& vFld
            ) const;

            //- Select candidate cells for refinement
            virtual void selectRefineCandidates
            (
                PackedBoolList& candidateCell,
                const scalar lowerRefineLevel,
                const scalar upperRefineLevel,
                const scalar maxRefinement,
                const scalarField& vFld,
                const labelList& cells
            ) const;

            //- Select candidate cells for refinement
            virtual scalar selectRefineCandidates
            (
                PackedBoolList& candidateCell,
                const dictionary& refineDict
            ) const;

            //- Subset candidate cells for refinement
            virtual labelList selectRefineCells
            (
                const label maxCells,
                const label maxRefinement,
                const PackedBoolList& candidateCell
            ) const;

            //- Get point max of connected cell
            void selectUnrefineCandidates
            (
                boolList& unrefineCandidates,
                const volScalarField&,
                const scalar unrefineLevel
            ) const;

            //- Get point max of connected cell
            void selectUnrefineCandidates
            (
                boolList& unrefineCandidates,
                const volScalarField&,
                const cellZone&,
                const scalar unrefineLevel
            ) const;

            //- Get point max of connected cell
            void selectUnrefineCandidates
            (
                boolList& unrefineCandidates,
                const dictionary& refineDict
            ) const;

            //- Select points that can be unrefined.
            virtual labelList selectUnrefinePoints
            (
                const PackedBoolList& markedCell,
                const boolList& unrefineCandidates
            ) const;

            //- Extend markedCell with cell-face-cell.
            void extendMarkedCells(PackedBoolList& markedCell) const;

            //- Check all cells have 8 anchor points
            void checkEightAnchorPoints
            (
                PackedBoolList& protectedCell,
                label& nProtected
            ) const;


public:

    //- Runtime type information
    TypeName("dynamicRefineFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit dynamicRefineFvMesh(const IOobject& io);

        //- Disallow default bitwise copy construction
        dynamicRefineFvMesh(const dynamicRefineFvMesh&) = delete;


    //- Destructor
    virtual ~dynamicRefineFvMesh();


    // Member Functions

        //- Direct access to the refinement engine
        const hexRef8& meshCutter() const
        {
            return meshCutter_;
        }

        //- Cells which should not be refined/unrefined
        const PackedBoolList& protectedCell() const
        {
            return protectedCells_;
        }

        //- Cells which should not be refined/unrefined
        PackedBoolList& protectedCell()
        {
            return protectedCells_;
        }

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();


    // Writing

        //- Write using given format, version and compression
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp,
            const bool write = true
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dynamicRefineFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
