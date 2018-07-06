/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::tabulatedWallFunctions::tabulatedWallFunction

Description
    Base class for models that generate tabulated wall function data.

SourceFiles
    tabulatedWallFunction.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedWallFunction_H
#define tabulatedWallFunction_H

#include "dictionary.H"
#include "polyMesh.H"
#include "runTimeSelectionTables.H"
#include "Switch.H"
#include "uniformInterpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace tabulatedWallFunctions
{

/*---------------------------------------------------------------------------*\
                   Class tabulatedWallFunction Declaration
\*---------------------------------------------------------------------------*/

class tabulatedWallFunction
{
protected:

    // Proteced data

        //- Main dictionary
        const dictionary dict_;

        //- Reference to the mesh database
        const polyMesh& mesh_;

        //- Coefficients dictionary
        const dictionary coeffDict_;

        //- Name of inverted table
        word invertedTableName_;

        //- Inverted table
        uniformInterpolationTable<scalar> invertedTable_;


public:

    //- Run-time type information
    TypeName("tabulatedWallFunction");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        tabulatedWallFunction,
        dictionary,
        (
            const dictionary& dict,
            const polyMesh& mesh
        ),
        (dict, mesh)
    );

    //- Constructor
    tabulatedWallFunction
    (
        const dictionary& dict,
        const polyMesh& mesh,
        const word& name
    );

    //- Destructor
    virtual ~tabulatedWallFunction();


    //- Selector
    static autoPtr<tabulatedWallFunction> New
    (
        const dictionary& dict,
        const polyMesh& mesh
    );


    // Member Functions

        // Access

            //- Return the inverted table name
            inline const word& invertedTableName() const;

            //- Return the inverted table
            inline const uniformInterpolationTable<scalar>&
                invertedTable() const;


        // I-O

            //- Write
            virtual void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace tabulatedWallFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabulatedWallFunctionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
