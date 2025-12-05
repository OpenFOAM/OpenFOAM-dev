/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2025 OpenFOAM Foundation
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
    Foam::pairGAMGProcAgglomeration

Description
    Pair-wise processor agglomerator.

    Neighbouring pairs of processors are agglomerated when the number of
    agglomerated cells per processor in the level is below \c
    minCellsPerProcessor and the pairs chosen based on the largest number of
    connecting processor faces.

Usage
    Example of GAMG solver settings with pair-wise processor agglomeration
    \verbatim
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;

        processorAgglomeration
        {
            agglomerator         pair;
            minCellsPerProcessor 300;
        }

        tolerance       1e-8;
        relTol          0;
    }
    \endverbatim

SourceFiles
    pairGAMGProcAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef pairGAMGProcAgglomeration_H
#define pairGAMGProcAgglomeration_H

#include "GAMGProcAgglomeration.H"
#include "DynamicList.H"
#include "labelField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class lduPrimitiveMesh;

/*---------------------------------------------------------------------------*\
               Class pairGAMGProcAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class pairGAMGProcAgglomeration
:
    public GAMGProcAgglomeration
{
    // Private Data

        //- When to processor agglomerate
        const label minCellsPerProcessor_;

        //- Allocated communicators
        DynamicList<label> comms_;

    // Private Member Functions

        //- Return (on master) all single-cell meshes collected. single-cell
        //  meshes are just one cell with all proc faces intact.
        autoPtr<lduPrimitiveMesh> singleCellMesh
        (
            const label singleCellMeshComm,
            const lduMesh& mesh,
            scalarField& faceWeights
        ) const;

        //- Construct processor agglomeration: for every processor the
        //  coarse processor-cluster it agglomerates onto
        tmp<labelField> processorAgglomeration(const lduMesh&) const;

        //- Do we need to agglomerate across processors?
        bool doProcessorAgglomeration(const lduMesh&) const;


public:

    //- Runtime type information
    TypeName("pair");


    // Constructors

        //- Construct given agglomerator and controls
        pairGAMGProcAgglomeration
        (
            GAMGAgglomeration& agglom,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        pairGAMGProcAgglomeration
        (
            const pairGAMGProcAgglomeration&
        ) = delete;


    //- Destructor
    virtual ~pairGAMGProcAgglomeration();


    // Member Functions

       //- Modify agglomeration. Return true if modified
        virtual bool agglomerate();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pairGAMGProcAgglomeration&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
