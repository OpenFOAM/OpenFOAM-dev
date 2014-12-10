/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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
    Foam::virtualMassModel

Description

SourceFiles
    virtualMassModel.C
    newVirtualMassModel.C

\*---------------------------------------------------------------------------*/

#ifndef virtualMassModel_H
#define virtualMassModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phasePair;

/*---------------------------------------------------------------------------*\
                         Class virtualMassModel Declaration
\*---------------------------------------------------------------------------*/

class virtualMassModel
:
    public regIOobject
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;


public:

    //- Runtime type information
    TypeName("virtualMassModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            virtualMassModel,
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

        //- Construct from a dictionary and a phase pair
        virtualMassModel
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~virtualMassModel();


    // Selectors

        static autoPtr<virtualMassModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Virtual mass coefficient
        virtual tmp<volScalarField> Cvm() const = 0;

        //- The virtual mass function K used in the momentum equation
        //    ddt(alpha1*rho1*U1) + ... = ... K*(DU1_Dt - DU2_Dt)
        //    ddt(alpha2*rho2*U2) + ... = ... K*(DU1_Dt - DU2_Dt)
        virtual tmp<volScalarField> K() const;

        // Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
