/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#ifndef ChemistryCombustion_H
#define ChemistryCombustion_H

#include "autoPtr.H"
#include "CombustionModel.H"
#include "BasicChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    class ChemistryCombustion Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class ChemistryCombustion
:
    public CombustionModel<ReactionThermo>
{
protected:

    // Protected data

        //- Pointer to chemistry model
        autoPtr<BasicChemistryModel<ReactionThermo>> chemistryPtr_;


public:

    // Constructors

        //- Construct from components and thermo
        ChemistryCombustion
        (
            const word& modelType,
            ReactionThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );


    //- Destructor
    virtual ~ChemistryCombustion();


    // Member Functions

        //- Return access to the thermo package
        virtual ReactionThermo& thermo();

        //- Return const access to the thermo package
        virtual const ReactionThermo& thermo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ChemistryCombustion.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
