/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::dragModels::timeScaleFiltered

Description
    A time scale filtering wrapper around an underlying drag model intended for
    simulations involving dispersed phase with a very small diameter.
    The model limits the drag coefficient based on the relaxation time given.

SourceFiles
    timeScaleFiltered.C

\*---------------------------------------------------------------------------*/

#ifndef timeScaleFilteredDrag_H
#define timeScaleFilteredDrag_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{

/*---------------------------------------------------------------------------*\
                    Class timeScaleFilteredDrag Declaration
\*---------------------------------------------------------------------------*/

class timeScaleFilteredDrag
:
    public dragModel
{
    // Private Data

        //- Pointer to the underlying drag model
        autoPtr<dragModel> dragModel_;

        //- The minimum relaxation time
        const dimensionedScalar minRelaxTime_;


public:

    //- Runtime type information
    TypeName("timeScaleFiltered");


    // Constructors

        //- Construct from a dictionary and a phase pair
        timeScaleFilteredDrag
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~timeScaleFilteredDrag();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
