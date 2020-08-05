/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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
    Foam::ChemistryCombustion

Description
    Chemistry model wrapper for combustion models

SourceFiles
    ChemistryCombustion.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryCombustion_H
#define chemistryCombustion_H

#include "autoPtr.H"
#include "combustionModel.H"
#include "basicChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                    class chemistryCombustion Declaration
\*---------------------------------------------------------------------------*/

class chemistryCombustion
:
    public combustionModel
{
protected:

    // Protected data

        //- Pointer to chemistry model
        autoPtr<basicChemistryModel> chemistryPtr_;


public:

    // Constructors

        //- Construct from components and thermo
        chemistryCombustion
        (
            const word& modelType,
            const fluidReactionThermo& thermo,
            const compressibleMomentumTransportModel& turb,
            const word& combustionProperties
        );


    //- Destructor
    virtual ~chemistryCombustion();


    // Member Functions

        //- Return const access to the chemistry
        const basicChemistryModel& chemistry() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
