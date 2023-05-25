/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::heatTransferModel

Description
    Model for heat transfer between phases

SourceFiles
    heatTransferModel.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferModel_H
#define heatTransferModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "BlendedInterfacialModel.H"
#include "SidedInterfacialModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class heatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class heatTransferModel
:
    public regIOobject
{
protected:

    // Protected data

        //- Residual phase fraction
        const dimensionedScalar residualAlpha_;


public:

    //- Runtime type information
    TypeName("heatTransferModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            heatTransferModel,
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

        //- Does this model require correcting on fixed flux boundaries?
        static const bool correctFixedFluxBCs = false;


    // Constructors

        //- Construct from a dictionary and an interface
        heatTransferModel
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~heatTransferModel();


    // Selectors

        static autoPtr<heatTransferModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool outer=true,
            const bool registerObject=true
        );


    // Member Functions

        //- The heat transfer function K used in the enthalpy equation
        //    ddt(alpha1*rho1*ha) + ... = ... K*(Ta - Tb)
        //    ddt(alpha2*rho2*hb) + ... = ... K*(Tb - Ta)
        tmp<volScalarField> K() const;

        //- The heat transfer function K used in the enthalpy equation
        //    ddt(alpha1*rho1*ha) + ... = ... K*(Ta - Tb)
        //    ddt(alpha2*rho2*hb) + ... = ... K*(Tb - Ta)
        //  with a specified residual volume fraction
        virtual tmp<volScalarField> K(const scalar residualAlpha) const = 0;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


/*---------------------------------------------------------------------------*\
                   Class blendedHeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class blendedHeatTransferModel
:
    public BlendedInterfacialModel<heatTransferModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using
            BlendedInterfacialModel<heatTransferModel>::
            BlendedInterfacialModel;


    // Selectors

        static autoPtr<blendedHeatTransferModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    // Member Functions

        //- Return the heatTransfer coefficient K
        tmp<volScalarField> K() const;

        //- Return the heatTransfer coefficient K
        tmp<volScalarField> K(const scalar residualAlpha) const;
};


/*---------------------------------------------------------------------------*\
                 Class sidedBlendedHeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class sidedBlendedHeatTransferModel
:
    public SidedInterfacialModel<blendedHeatTransferModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using
            SidedInterfacialModel<blendedHeatTransferModel>::
            SidedInterfacialModel;


    // Selectors

        static autoPtr<sidedBlendedHeatTransferModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
