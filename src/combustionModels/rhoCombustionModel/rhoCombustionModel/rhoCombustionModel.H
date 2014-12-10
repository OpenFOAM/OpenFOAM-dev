/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::rhoCombustionModel

Description
    Combustion models for rho-based thermodynamics

SourceFiles
    rhoCombustionModelI.H
    rhoCombustionModel.C
    rhoCombustionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef rhoCombustionModel_H
#define rhoCombustionModel_H

#include "combustionModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                     class rhoCombustionModel Declaration
\*---------------------------------------------------------------------------*/

class rhoCombustionModel
:
    public combustionModel
{

    // Private Member Functions

        //- Construct as copy (not implemented)
        rhoCombustionModel(const rhoCombustionModel&);

        //- Disallow default bitwise assignment
        void operator=(const rhoCombustionModel&);


public:

    //- Runtime type information
    TypeName("rhoCombustionModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoCombustionModel,
        dictionary,
        (
            const word& modelType,
            const fvMesh& mesh
        ),
        (modelType, mesh)
    );


    // Constructors


        //- Construct from components and thermo
        rhoCombustionModel(const word& modelType, const fvMesh& mesh);



    //- Selector
    static autoPtr<rhoCombustionModel> New(const fvMesh& mesh);


    //- Destructor
    virtual ~rhoCombustionModel();


    // Member Functions

        // Access functions

            //- Access combustion dict
             inline const dictionary& coeff() const;


        //- Return access to the thermo package
        virtual rhoReactionThermo& thermo() = 0;

        //- Return const access to the thermo package
        virtual const rhoReactionThermo& thermo() const = 0;

        //- Return tmp of rho
        virtual tmp<volScalarField> rho() const = 0;


     // I-O

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
