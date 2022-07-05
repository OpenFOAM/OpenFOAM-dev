/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::fvMeshTopoChangers::refiner

Description
    Dynamic mesh refinement/unrefinement based on volScalarField values.

    Refinement can optionally be specified in a cellZone or in multiple
    regions, each controlled by a different volScalarField.

Usage
    Example of single field based refinement in all cells:
    \verbatim
    topoChanger
    {
        type            refiner;

        libs            ("libfvMeshTopoChangers.so");

        // How often to refine
        refineInterval  1;

        // Field to be refinement on
        field           alpha.water;

        // Refine field in between lower..upper
        lowerRefineLevel 0.001;
        upperRefineLevel 0.999;

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
    }
    \endverbatim

    Example of single field based refinement in two regions:
    \verbatim
    topoChanger
    {
        type            refiner;

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
    }
    \endverbatim


SourceFiles
    fvMeshTopoChangersRefiner.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshTopoChangersRefiner_H
#define fvMeshTopoChangersRefiner_H

#include "fvMeshTopoChanger.H"
#include "hexRef8.H"
#include "PackedBoolList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{

/*---------------------------------------------------------------------------*\
                 Class fvMeshTopoChangers::refiner Declaration
\*---------------------------------------------------------------------------*/

class refiner
:
    public fvMeshTopoChanger
{
    // Private Member Data

        //- Refinement control dictionary
        dictionary dict_;

        label refineInterval_;

        label maxCells_;

        label nBufferLayers_;

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

        //- The time index used for updating
        label timeIndex_;


    // Private Member Functions

        //- Count set/unset elements in packedlist.
        static label count(const PackedBoolList&, const unsigned int);

        //- Calculate cells that cannot be refined since would trigger
        //  refinement of protectedCells_ (since 2:1 refinement cascade)
        void calculateProtectedCells(PackedBoolList& unrefineableCells) const;

        //- Read the parameters from dictionary
        void readDict();

        //- Refine cells. Update mesh and fields.
        autoPtr<polyTopoChangeMap> refine(const labelList&);

        //- Unrefine cells. Gets passed in centre points of cells to combine.
        autoPtr<polyTopoChangeMap> unrefine(const labelList&);

        //- Find the U field name corresponding to Uf
        word Uname(const surfaceVectorField& Uf) const;

        void refineFluxes
        (
            const labelHashSet& masterFaces,
            const polyTopoChangeMap& map
        );

        void refineUfs
        (
            const labelHashSet& masterFaces,
            const polyTopoChangeMap& map
        );

        void unrefineFluxes
        (
            const Map<label>& faceToSplitPoint,
            const polyTopoChangeMap& map
        );

        void unrefineUfs
        (
            const Map<label>& faceToSplitPoint,
            const polyTopoChangeMap& map
        );


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

            //- Select points that can be unrefined.
            virtual labelList selectUnrefinePoints
            (
                const PackedBoolList& markedCell
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
    TypeName("refiner");


    // Constructors

        //- Construct from fvMesh and dictionary
        refiner(fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        refiner(const refiner&) = delete;


    //- Destructor
    virtual ~refiner();


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

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Writing

        //- Write using given format, version and compression
        virtual bool write(const bool write = true) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const refiner&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshTopoChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
