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
    Foam::regionModels::surfaceFilmModels::transferModelList

Description
    List container for film transfer models

SourceFiles
    transferModelList.C

\*---------------------------------------------------------------------------*/

#ifndef transferModelList_H
#define transferModelList_H

#include "PtrList.H"
#include "transferModel.H"
#include "filmSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class transferModelList Declaration
\*---------------------------------------------------------------------------*/

class transferModelList
:
    public PtrList<transferModel>,
    public filmSubModelBase
{
    // Private data

        //- List of mass transferred per patch
        scalarField massTransferred_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        transferModelList(const transferModelList&);

        //- Disallow default bitwise assignment
        void operator=(const transferModelList&);


public:

    // Constructors

        //- Construct null
        transferModelList(surfaceFilmRegionModel& film);

        //- Construct from type name, dictionary and surface film model
        transferModelList
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~transferModelList();


    // Member Functions

        //- Correct kinematic transfers
        virtual void correct
        (
            scalarField& availableMass,
            volScalarField& massToTransfer
        );

        //- Correct kinematic and thermodynamic transfers
        virtual void correct
        (
            scalarField& availableMass,
            volScalarField& massToTransfer,
            volScalarField& energyToTransfer
        );

        //- Provide some info
        virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
