/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2025 OpenFOAM Foundation
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
    Foam::populationBalance::shapeModel

Description
    Base class for modelling the shape of the particles belonging to a size
    class through alternative diameters, e.g. a collisional diameter, which can
    then be utilised in population balance submodels, e.g. for modelling
    fractal aggregation.

SourceFiles
    shapeModel.C

\*---------------------------------------------------------------------------*/

#ifndef shapeModel_H
#define shapeModel_H

#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class populationBalanceModel;

namespace populationBalance
{

/*---------------------------------------------------------------------------*\
                         Class shapeModel Declaration
\*---------------------------------------------------------------------------*/

class shapeModel
{
protected:

    // Protected Data

        //- Reference to the population balance model
        const populationBalanceModel& popBal_;


public:

    //- Runtime type information
    TypeName("shapeModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            shapeModel,
            dictionary,
            (
                const dictionary& dict,
                const populationBalanceModel& popBal
            ),
            (dict, popBal)
        );


    // Constructors

        //- Construct from a population balance model
        shapeModel(const populationBalanceModel& popBal);

        //- Disallow default bitwise copy construction
        shapeModel(const shapeModel&) = delete;


    // Selectors

        static autoPtr<shapeModel> New
        (
            const dictionary& dict,
            const populationBalanceModel& popBal
        );


    //- Destructor
    virtual ~shapeModel();


    // Member Functions

        // Access

            //- Return reference to the population balance model
            const populationBalanceModel& popBal() const;

            //- Return representative surface area of the group
            virtual tmp<volScalarField> a(const label i) const = 0;

            //- Return representative diameter of the group
            virtual tmp<volScalarField> d(const label i) const = 0;


        // Edit

            //- Solve equations
            virtual void solve();

            //- Correct the shape
            virtual void correct();

            //- Add coalescence contribution to secondary property source
            virtual void addCoalescence
            (
                const volScalarField::Internal& Su,
                const label i,
                const label j,
                const label k
            );

            //- Add breakup contribution to secondary property source
            virtual void addBreakup
            (
                const volScalarField::Internal& Su,
                const label i,
                const label j
            );

            //- Reset secondary property source
            virtual void reset();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const shapeModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalance
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
