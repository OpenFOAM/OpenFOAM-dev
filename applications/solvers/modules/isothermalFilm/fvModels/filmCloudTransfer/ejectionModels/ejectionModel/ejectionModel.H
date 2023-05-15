/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::filmEjectionModels::ejectionModel

Description
    Abstract base class for film to cloud ejection transfer models

SourceFiles
    ejectionModel.C
    ejectionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef ejectionModel_H
#define ejectionModel_H

#include "isothermalFilm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ejectionModel Declaration
\*---------------------------------------------------------------------------*/

class ejectionModel
{

protected:

    // Protected Member Data

        //- Reference to the film
        const solvers::isothermalFilm& film_;

        //- Cached ejection rate (transfer to cloud)
        volScalarField::Internal rate_;

        //- Cached diameter of ejected droplets
        volScalarField::Internal diameter_;


public:

    //- Runtime type information
    TypeName("ejectionModel");


    // Declare run-time constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             ejectionModel,
             dictionary,
             (
                 const dictionary& dict,
                 const solvers::isothermalFilm& film
             ),
             (dict, film)
         );


    // Constructors

        //- Construct from dictionary and film model
        ejectionModel
        (
            const dictionary& dict,
            const solvers::isothermalFilm& film
        );

        //- Disallow default bitwise copy construction
        ejectionModel(const ejectionModel&) = delete;


    // Selector
    static autoPtr<ejectionModel> New
    (
        const dictionary& dict,
        const solvers::isothermalFilm& film
    );


    //- Destructor
    virtual ~ejectionModel();


    // Member Functions

        //- Return ejection rate (transfer to cloud)
        const volScalarField::Internal& rate() const
        {
            return rate_;
        }

        //- Return diameter of ejected droplets
        const volScalarField::Internal& diameter() const
        {
            return diameter_;
        }

        //- Correct
        virtual void correct() = 0;


        // Mesh changes

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ejectionModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
