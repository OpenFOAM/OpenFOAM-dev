/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::dragModel

Description

SourceFiles
    dragModel.C
    newDragModel.C

\*---------------------------------------------------------------------------*/

#ifndef dragModel_H
#define dragModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phasePair;
class swarmCorrection;

/*---------------------------------------------------------------------------*\
                           Class dragModel Declaration
\*---------------------------------------------------------------------------*/

class dragModel
:
    public regIOobject
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;

        //- Swarm correction
        autoPtr<swarmCorrection> swarmCorrection_;

        //- Residual phase fraction
        const dimensionedScalar residualAlpha_;


public:

    //- Runtime type information
    TypeName("dragModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            dragModel,
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
        static const dimensionSet dimK;


    // Constructors

        // Construct without residual constants
        dragModel
        (
            const phasePair& pair,
            const bool registerObject
        );

        // Construct with residual constants
        dragModel
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~dragModel();


    // Selectors

        static autoPtr<dragModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const = 0;

        //- The drag function K used in the momentum equation
        //    ddt(alpha1*rho1*U1) + ... = ... K*(U1-U2)
        //    ddt(alpha2*rho2*U2) + ... = ... K*(U2-U1)
        virtual tmp<volScalarField> K() const;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
