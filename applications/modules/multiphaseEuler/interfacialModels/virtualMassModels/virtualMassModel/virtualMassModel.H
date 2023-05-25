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
    Foam::virtualMassModel

Description
    Model for virtual mass between phases

SourceFiles
    virtualMassModel.C
    virtualMassModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef virtualMassModel_H
#define virtualMassModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "BlendedInterfacialModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class virtualMassModel Declaration
\*---------------------------------------------------------------------------*/

class virtualMassModel
:
    public regIOobject
{
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
                const phaseInterface& interface,
                const bool registerObject
            ),
            (dict, interface, registerObject)
        );


    // Static Data Members

        //- Coefficient dimensions
        static const dimensionSet dimK;

        //- This model should not be set to 0 on fixed flux boundaries
        static const bool correctFixedFluxBCs = false;


    // Constructors

        //- Construct from a dictionary and an interface
        virtualMassModel
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~virtualMassModel();


    // Selectors

        static autoPtr<virtualMassModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool outer=true,
            const bool registerObject=true
        );


    // Member Functions

        //- Return the virtual mass coefficient K
        //  used in the momentum equation
        //    ddt(alpha1*rho1*U1) + ... = ... K*(DU1_Dt - DU2_Dt)
        //    ddt(alpha2*rho2*U2) + ... = ... K*(DU1_Dt - DU2_Dt)
        virtual tmp<volScalarField> K() const = 0;

        // Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


/*---------------------------------------------------------------------------*\
                      Class blendedVirtualMassModel Declaration
\*---------------------------------------------------------------------------*/

class blendedVirtualMassModel
:
    public BlendedInterfacialModel<virtualMassModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using
            BlendedInterfacialModel<virtualMassModel>
          ::BlendedInterfacialModel;


    // Selectors

        static autoPtr<blendedVirtualMassModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    // Member Functions

        //- Return the virtual mass coefficient K
        tmp<volScalarField> K() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
