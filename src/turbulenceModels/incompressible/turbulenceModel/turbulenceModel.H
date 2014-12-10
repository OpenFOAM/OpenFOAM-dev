/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Namespace
    Foam::incompressible::turbulenceModels

Description
    Namespace for incompressible turbulence turbulence models.

Class
    Foam::incompressible::turbulenceModel

Description
    Abstract base class for incompressible turbulence models
    (RAS, LES and laminar).

SourceFiles
    turbulenceModel.C
    newTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef turbulenceModel_H
#define turbulenceModel_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "incompressible/transportModel/transportModel.H"
#include "nearWallDist.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

namespace incompressible
{

/*---------------------------------------------------------------------------*\
                           Class turbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class turbulenceModel
:
    public regIOobject
{

protected:

    // Protected data

        const Time& runTime_;
        const fvMesh& mesh_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        transportModel& transportModel_;

        //- Near wall distance boundary field
        nearWallDist y_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        turbulenceModel(const turbulenceModel&);

        //- Disallow default bitwise assignment
        void operator=(const turbulenceModel&);


public:

    //- Runtime type information
    TypeName("turbulenceModel");


    // Declare run-time New selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            turbulenceModel,
            turbulenceModel,
            (
                const volVectorField& U,
                const surfaceScalarField& phi,
                transportModel& transport,
                const word& turbulenceModelName
            ),
            (U, phi, transport, turbulenceModelName)
        );


    // Constructors

        //- Construct from components
        turbulenceModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = typeName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<turbulenceModel> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = typeName
        );


    //- Destructor
    virtual ~turbulenceModel()
    {}


    // Member Functions

        //- Const access to the coefficients dictionary
        virtual const dictionary& coeffDict() const = 0;

        //- Helper function to return the nam eof the turbulence G field
        inline word GName() const
        {
            return word(type() + ":G");
        }

        //- Access function to velocity field
        inline const volVectorField& U() const
        {
            return U_;
        }

        //- Access function to flux field
        inline const surfaceScalarField& phi() const
        {
            return phi_;
        }

        //- Access function to incompressible transport model
        inline transportModel& transport() const
        {
            return transportModel_;
        }

        //- Return the near wall distances
        const nearWallDist& y() const
        {
            return y_;
        }

        //- Return the laminar viscosity
        inline tmp<volScalarField> nu() const
        {
            return transportModel_.nu();
        }

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const = 0;

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const = 0;

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const = 0;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const = 0;

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const = 0;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const = 0;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const = 0;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const = 0;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct() = 0;

        //- Read LESProperties or RASProperties dictionary
        virtual bool read() = 0;

        //- Default dummy write function
        virtual bool writeData(Ostream&) const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
