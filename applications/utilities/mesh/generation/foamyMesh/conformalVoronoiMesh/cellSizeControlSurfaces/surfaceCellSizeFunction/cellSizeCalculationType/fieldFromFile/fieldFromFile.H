/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Foam::fieldFromFile

Description

SourceFiles
    fieldFromFile.C

\*---------------------------------------------------------------------------*/

#ifndef fieldFromFile_H
#define fieldFromFile_H

#include "cellSizeCalculationType.H"
#include "triSurfaceFields.H"
#include "PrimitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class triSurfaceMesh;

/*---------------------------------------------------------------------------*\
                       Class fieldFromFile Declaration
\*---------------------------------------------------------------------------*/

class fieldFromFile
:
    public cellSizeCalculationType
{

private:

    // Private data

        //- Dictionary of coefficients for automatic cell sizing
        const dictionary& coeffsDict_;

        //- Name of the triSurfaceScalarField file to load in. Must be in
        //  constant/triSurface
        const word fileName_;

        //- Multiply all the point sizes by this value
        const scalar cellSizeMultipleCoeff_;


public:

    //- Runtime type information
    TypeName("fieldFromFile");


    // Constructors

        //- Construct from components
        fieldFromFile
        (
            const dictionary& cellSizeCalcTypeDict,
            const triSurfaceMesh& surface,
            const scalar& defaultCellSize
        );


    //- Destructor
    virtual ~fieldFromFile()
    {}


    // Member Functions

        //- Load the cell size field
        virtual tmp<triSurfacePointScalarField> load();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
