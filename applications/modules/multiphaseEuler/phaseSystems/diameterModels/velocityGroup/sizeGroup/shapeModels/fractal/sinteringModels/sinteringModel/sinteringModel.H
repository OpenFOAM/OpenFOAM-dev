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
    Foam::diameterModels::shapeModels::sinteringModel

Description
    Abstract base class for modelling sintering of primary particles in fractal
    aggregates.

See also
    Foam::diameterModels::shapeModels::fractal

SourceFiles
    sinteringModel.C

\*---------------------------------------------------------------------------*/

#ifndef sinteringModel_H
#define sinteringModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "fractal.H"

namespace Foam
{
namespace diameterModels
{
namespace shapeModels
{

/*---------------------------------------------------------------------------*\
                       Class sinteringModel Declaration
\*---------------------------------------------------------------------------*/

class sinteringModel
{
protected:

    // Protected Data

        //- Reference to fractal shapeModel class
        const fractal& fractal_;


public:

    //- Runtime type information
    TypeName("sinteringModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            sinteringModel,
            dictionary,
            (
                const dictionary& dict,
                const fractal& fractalShape
            ),
            (dict, fractalShape)
        );


    // Constructors

        //- Construct from a dictionary and a fractal shape model
        sinteringModel
        (
            const dictionary& dict,
            const fractal& fractalShape
        );


    //- Destructor
    virtual ~sinteringModel();


    // Selectors

        static autoPtr<sinteringModel> New
        (
            const dictionary& dict,
            const fractal& fractalShape
        );


    // Member Functions

        //- Sintering source term
        virtual tmp<fvScalarMatrix> R() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace shapeModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
