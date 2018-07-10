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

Class
    Foam::cellSizeCalculationType

Description
    Abstract base class for specifying target cell sizes

SourceFiles
    cellSizeCalculationType.C

\*---------------------------------------------------------------------------*/

#ifndef cellSizeCalculationType_H
#define cellSizeCalculationType_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "triSurfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class triSurfaceMesh;

/*---------------------------------------------------------------------------*\
                     Class cellSizeCalculationType Declaration
\*---------------------------------------------------------------------------*/

class cellSizeCalculationType
{

protected:

        const dictionary& cellSizeCalculationTypeDict_;

        //- Reference to the triSurfaceMesh
        const triSurfaceMesh& surface_;

        const scalar& defaultCellSize_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        cellSizeCalculationType(const cellSizeCalculationType&);

        //- Disallow default bitwise assignment
        void operator=(const cellSizeCalculationType&);


public:

    //- Runtime type information
    TypeName("cellSizeCalculationType");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            cellSizeCalculationType,
            dictionary,
            (
                const dictionary& cellSizeCalculationTypeDict,
                const triSurfaceMesh& surface,
                const scalar& defaultCellSize
            ),
            (cellSizeCalculationTypeDict, surface, defaultCellSize)
        );


    // Constructors

        //- Construct from components
        cellSizeCalculationType
        (
            const word& type,
            const dictionary& cellSizeCalculationTypeDict,
            const triSurfaceMesh& surface,
            const scalar& defaultCellSize
        );


    // Selectors

        //- Return a reference to the selected cellSizeCalculationType
        static autoPtr<cellSizeCalculationType> New
        (
            const dictionary& cellSizeCalculationTypeDict,
            const triSurfaceMesh& surface,
            const scalar& defaultCellSize
        );


    //- Destructor
    virtual ~cellSizeCalculationType();


    // Member Functions

        //- Load the cell size
        virtual tmp<triSurfacePointScalarField> load() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
