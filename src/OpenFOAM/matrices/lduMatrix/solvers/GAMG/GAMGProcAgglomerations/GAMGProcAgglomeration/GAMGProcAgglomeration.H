/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::GAMGProcAgglomeration

Description
    Processor agglomeration of GAMGAgglomerations.

SourceFiles
    GAMGProcAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef GAMGProcAgglomeration_H
#define GAMGProcAgglomeration_H

#include "runTimeSelectionTables.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class GAMGAgglomeration;
class lduMesh;

/*---------------------------------------------------------------------------*\
                    Class GAMGProcAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class GAMGProcAgglomeration
{

protected:

    // Protected data

        //- Reference to agglomeration
        GAMGAgglomeration& agglom_;

    // Protected Member Functions

        //- Debug: write agglomeration info
        void printStats(Ostream& os, GAMGAgglomeration& agglom) const;

        //- Agglomerate a level. Return true if anything has changed
        bool agglomerate
        (
            const label fineLevelIndex,
            const labelList& procAgglomMap,
            const labelList& masterProcs,
            const List<int>& agglomProcIDs,
            const label procAgglomComm
        );

        //- Debug: calculate global cell-cells
        static labelListList globalCellCells(const lduMesh&);

private:

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        GAMGProcAgglomeration(const GAMGProcAgglomeration&);

        //- Disallow default bitwise assignment
        void operator=(const GAMGProcAgglomeration&);


public:

    //- Runtime type information
    TypeName("GAMGProcAgglomeration");


    // Declare run-time constructor selection tables

        //- Runtime selection table for pure geometric agglomerators
        declareRunTimeSelectionTable
        (
            autoPtr,
            GAMGProcAgglomeration,
            GAMGAgglomeration,
            (
                GAMGAgglomeration& agglom,
                const dictionary& controlDict
            ),
            (
                agglom,
                controlDict
            )
        );


    // Constructors

        //- Construct given agglomerator and controls
        GAMGProcAgglomeration
        (
            GAMGAgglomeration& agglom,
            const dictionary& controlDict
        );


    // Selectors

        //- Return the selected agglomerator
        static autoPtr<GAMGProcAgglomeration> New
        (
            const word& type,
            GAMGAgglomeration& agglom,
            const dictionary& controlDict
        );


    //- Destructor
    virtual ~GAMGProcAgglomeration();


    // Member Functions

        //- Modify agglomeration. Return true if modified
        virtual bool agglomerate() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
