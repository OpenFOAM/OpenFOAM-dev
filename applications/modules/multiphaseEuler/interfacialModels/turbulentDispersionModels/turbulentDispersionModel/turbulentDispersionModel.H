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
    Foam::turbulentDispersionModel

Description
    Model for turbulent dispersion between two phases

SourceFiles
    turbulentDispersionModel.C
    turbulentDispersionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentDispersionModel_H
#define turbulentDispersionModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "BlendedInterfacialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class turbulentDispersionModel Declaration
\*---------------------------------------------------------------------------*/

class turbulentDispersionModel
{
public:

    //- Runtime type information
    TypeName("turbulentDispersionModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            turbulentDispersionModel,
            dictionary,
            (
                const dictionary& dict,
                const phaseInterface& interface
            ),
            (dict, interface)
        );

    // Static Data Members

        //- Diffusivity dimensions
        static const dimensionSet dimD;

        //- Does this model require correcting on fixed flux boundaries?
        static const bool correctFixedFluxBCs = true;


    // Constructors

        //- Construct from a dictionary and an interface
        turbulentDispersionModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~turbulentDispersionModel();


    // Selectors

        static autoPtr<turbulentDispersionModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool outer=true
        );


    // Member Functions

        //- Return the turbulent diffusivity multiplying the gradient of the
        //  phase-fraction
        virtual tmp<volScalarField> D() const = 0;
};


/*---------------------------------------------------------------------------*\
                Class blendedTurbulentDispersionModel Declaration
\*---------------------------------------------------------------------------*/

class blendedTurbulentDispersionModel
:
    public BlendedInterfacialModel<turbulentDispersionModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using
            BlendedInterfacialModel<turbulentDispersionModel>::
            BlendedInterfacialModel;


    // Selectors

        static autoPtr<blendedTurbulentDispersionModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    // Member Functions

        //- Return the turbulent diffusivity multiplying the gradient of the
        //  phase-fraction
        tmp<volScalarField> D() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
