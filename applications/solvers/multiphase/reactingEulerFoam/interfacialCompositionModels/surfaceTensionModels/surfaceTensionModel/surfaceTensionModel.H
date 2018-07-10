/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::surfaceTensionModel

Description

SourceFiles
    surfaceTensionModel.C
    newAspectRatioModel.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceTensionModel_H
#define surfaceTensionModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phasePair;

/*---------------------------------------------------------------------------*\
                           Class surfaceTensionModel Declaration
\*---------------------------------------------------------------------------*/

class surfaceTensionModel
:
    public regIOobject
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;


public:

    //- Runtime type information
    TypeName("surfaceTensionModel");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        surfaceTensionModel,
        dictionary,
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        ),
        (dict, pair, registerObject)
    );


    // Static data members

        //- Coefficient dimensions
        static const dimensionSet dimSigma;


    // Constructors

        //- Construct from a dictionary and a phase pair
        surfaceTensionModel
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~surfaceTensionModel();


    // Selectors

        static autoPtr<surfaceTensionModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Aspect ratio
        virtual tmp<volScalarField> sigma() const = 0;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
