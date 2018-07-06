/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::transferModel

Description
    Base class for film transfer models, handling mass transfer between the
    film and the continuous phase.

SourceFiles
    transferModel.C
    transferModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef transferModel_H
#define transferModel_H

#include "filmSubModelBase.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class transferModel Declaration
\*---------------------------------------------------------------------------*/

class transferModel
:
    public filmSubModelBase
{
    // Private data

        //- Transferred mass
        scalar transferredMass_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        transferModel(const transferModel&);

        //- Disallow default bitwise assignment
        void operator=(const transferModel&);


protected:

    // Protected Member Functions

        //- Add to transferred mass
        void addToTransferredMass(const scalar dMass);

        //- Correct
        void correct();


public:

    //- Runtime type information
    TypeName("transferModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             transferModel,
             dictionary,
             (
                surfaceFilmRegionModel& film,
                const dictionary& dict
             ),
             (film, dict)
         );


    // Constructors

        //- Construct for film
        transferModel(surfaceFilmRegionModel& film);

        //- Construct from type name, dictionary and surface film model
        transferModel
        (
            const word& modelType,
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected injection model
        static autoPtr<transferModel> New
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict,
            const word& modelType
        );


    //- Destructor
    virtual ~transferModel();


    // Member Functions

        //- Correct kinematic transfers
        virtual void correct
        (
            scalarField& availableMass,
            scalarField& massToTransfer
        ) = 0;

        //- Correct kinematic and thermodynamic transfers
        virtual void correct
        (
            scalarField& availableMass,
            scalarField& massToTransfer,
            scalarField& energyToTransfer
        );

        //- Return the total mass transferred
        virtual scalar transferredMassTotal() const;

        //- Accumulate the total mass transferred for the patches into the
        //  scalarField provided
        virtual void patchTransferredMassTotals(scalarField& patchMasses) const
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
