/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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
    Foam::diameterModels::shapeModel

Description
    Base class for modelling the shape of the particles belonging to a size
    class through alternative diameters, e.g. a collisional diameter, which can
    then be utilised in population balance submodels, e.g. for modelling
    fractal aggregation.

See also
    Foam::diameterModels::sizeGroup

SourceFiles
    shapeModel.C

\*---------------------------------------------------------------------------*/

#ifndef shapeModel_H
#define shapeModel_H

#include "dictionary.H"
#include "sizeGroup.H"
#include "driftModel.H"
#include "nucleationModel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                         Class shapeModel Declaration
\*---------------------------------------------------------------------------*/

class shapeModel
{

protected:

    // Protected Data

        //- Reference to sizeGroup
        const sizeGroup& sizeGroup_;


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
                const sizeGroup& group
            ),
            (dict, group)
        );


    // Constructors

        //- Construct from dictionary and sizeGroup
        shapeModel
        (
            const dictionary& dict,
            const sizeGroup& group
        );

        //- Disallow default bitwise copy construction
        shapeModel(const shapeModel&) = delete;


    // Selectors

        static autoPtr<shapeModel> New
        (
            const dictionary& dict,
            const sizeGroup& group
        );


    //- Destructor
    virtual ~shapeModel();


    // Member Functions

        // Access

            //- Return reference to size group
            const sizeGroup& SizeGroup() const;

            //- Return representative surface area of the sizeGroup
            virtual const tmp<volScalarField> a() const = 0;

            //- Return representative diameter of the sizeGroup
            virtual const tmp<volScalarField> d() const = 0;

        // Edit

            //- Correct physical diameter
            virtual void correct();

            //- Add coalescence contribution to secondary property source
            virtual void addCoalescence
            (
                const volScalarField& Su,
                const sizeGroup& fj,
                const sizeGroup& fk
            );

            //- Add breakup contribution to secondary property source
            virtual void addBreakup
            (
                const volScalarField& Su,
                const sizeGroup& fj
            );

            //- Add drift contribution to secondary property source
            virtual void addDrift
            (
                const volScalarField& Su,
                const sizeGroup& fu,
                const driftModel& model
            );

            //- Add nucleation contribution to secondary property source
            virtual void addNucleation
            (
                const volScalarField& Su,
                const sizeGroup& fi,
                const nucleationModel& model
            );

            //- Reset secondary property source
            virtual void reset();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const shapeModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
