/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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
    Foam::aspectRatioModel

Description
    Model for deviations in the shape of the dispersed phase from spherical.
    Just a sub-model modifier, typically for the drag model. Not a proper part
    of the diameter/shape modelling in the phase models.

SourceFiles
    aspectRatioModel.C
    aspectRatioModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef aspectRatioModel_H
#define aspectRatioModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "dispersedPhaseInterface.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class aspectRatioModel Declaration
\*---------------------------------------------------------------------------*/

class aspectRatioModel
{
protected:

    // Protected data

        //- Interface
        const dispersedPhaseInterface interface_;


public:

    //- Runtime type information
    TypeName("aspectRatioModel");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        aspectRatioModel,
        dictionary,
        (
            const dictionary& dict,
            const phaseInterface& interface
        ),
        (dict, interface)
    );


    // Constructors

        //- Construct from a dictionary and an interface
        aspectRatioModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~aspectRatioModel();


    // Selectors

        static autoPtr<aspectRatioModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    // Member Functions

        //- Aspect ratio
        virtual tmp<volScalarField> E() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
