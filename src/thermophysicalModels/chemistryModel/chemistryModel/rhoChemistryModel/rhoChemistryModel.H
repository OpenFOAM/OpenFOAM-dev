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
    Foam::rhoChemistryModel

Description
    Chemistry model for density-based thermodynamics

SourceFiles
    rhoChemistryModelI.H
    rhoChemistryModel.C
    newChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef rhoChemistryModel_H
#define rhoChemistryModel_H

#include "basicChemistryModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class rhoChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class rhoChemistryModel
:
    public basicChemistryModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        rhoChemistryModel(const rhoChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const rhoChemistryModel&);


protected:

    // Protected data

        //- Thermo package
        autoPtr<rhoReactionThermo> thermo_;


public:

    //- Runtime type information
    TypeName("rho");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoChemistryModel,
        fvMesh,
        (
            const fvMesh& mesh
        ),
        (mesh)
    );


    // Constructors

        //- Construct from mesh
        rhoChemistryModel(const fvMesh& mesh);


    //- Selector
    static autoPtr<rhoChemistryModel> New(const fvMesh& mesh);


    //- Destructor
    virtual ~rhoChemistryModel();


    // Member Functions

        //- Return access to the thermo package
        inline rhoReactionThermo& thermo();

        //- Return const access to the thermo package
        inline const rhoReactionThermo& thermo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rhoChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
