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
    Foam::psiChemistryModel

Description
    Chemistry model for compressibility-based thermodynamics

SourceFiles
    psiChemistryModelI.H
    psiChemistryModel.C
    newChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef psiChemistryModel_H
#define psiChemistryModel_H

#include "basicChemistryModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "psiReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class psiChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class psiChemistryModel
:
    public basicChemistryModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        psiChemistryModel(const psiChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const psiChemistryModel&);


protected:

    // Protected data

        //- Thermo package
        autoPtr<psiReactionThermo> thermo_;


public:

    //- Runtime type information
    TypeName("psi");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        psiChemistryModel,
        fvMesh,
        (const fvMesh& mesh),
        (mesh)
    );


    // Constructors

        //- Construct from mesh
        psiChemistryModel(const fvMesh& mesh);


    //- Selector
    static autoPtr<psiChemistryModel> New(const fvMesh& mesh);


    //- Destructor
    virtual ~psiChemistryModel();


    // Member Functions

        //- Return access to the thermo package
        inline psiReactionThermo& thermo();

        //- Return const access to the thermo package
        inline const psiReactionThermo& thermo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "psiChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
