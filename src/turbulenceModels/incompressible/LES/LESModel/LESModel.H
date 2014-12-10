/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::incompressible::LESModels

Group
    grpIcoLESTurbulence

Description
    Namespace for incompressible LES models.

Class
    Foam::incompressible::LESModel

Description
    Base class for all incompressible flow LES SGS models.

    This class defines the basic interface for an incompressible flow SGS
    model, and encapsulates data of value to all possible models.
    In particular this includes references to all the dependent fields
    (U, phi), the physical viscosity nu, and the LESProperties
    dictionary, which contains the model selection and model coefficients.

SourceFiles
    LESModel.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleLESModel_H
#define incompressibleLESModel_H

#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "incompressible/LES/incompressibleLESdelta/incompressibleLESdelta.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "incompressible/transportModel/transportModel.H"
#include "bound.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
                           Class LESModel Declaration
\*---------------------------------------------------------------------------*/

class LESModel
:
    public turbulenceModel,
    public IOdictionary
{

protected:

    // Protected data

        Switch printCoeffs_;
        dictionary coeffDict_;

        dimensionedScalar kMin_;

        autoPtr<Foam::LESdelta> delta_;


    // Protected Member Functions

        //- Print model coefficients
        virtual void printCoeffs();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        LESModel(const LESModel&);

        //- Disallow default bitwise assignment
        LESModel& operator=(const LESModel&);


public:

    //- Runtime type information
    TypeName("LESModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            LESModel,
            dictionary,
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
        LESModel
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName
        );


    // Selectors

        //- Return a reference to the selected LES model
        static autoPtr<LESModel> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName
        );


    //- Destructor
    virtual ~LESModel()
    {}


    // Member Functions

        // Access

            //- Const access to the coefficients dictionary,
            //  which provides info. about choice of models,
            //  and all related data (particularly model coefficients).
            virtual const dictionary& coeffDict() const
            {
                return coeffDict_;
            }

            //- Return the lower allowable limit for k (default: SMALL)
            const dimensionedScalar& kMin() const
            {
                return kMin_;
            }

            //- Allow kMin to be changed
            dimensionedScalar& kMin()
            {
                return kMin_;
            }

            //- Access function to filter width
            virtual const volScalarField& delta() const
            {
                return delta_();
            }


        //- Return the SGS viscosity.
        virtual tmp<volScalarField> nuSgs() const = 0;

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("nuEff", nuSgs() + nu())
            );
        }

        //- Return the sub-grid stress tensor.
        virtual tmp<volSymmTensorField> B() const = 0;


        // RAS compatibility functions for the turbulenceModel base class

            //- Return the turbulence viscosity
            virtual tmp<volScalarField> nut() const
            {
                return nuSgs();
            }

            //- Return the Reynolds stress tensor
            virtual tmp<volSymmTensorField> R() const
            {
                return B();
            }


        //- Correct Eddy-Viscosity and related properties.
        //  This calls correct(const tmp<volTensorField>& gradU) by supplying
        //  gradU calculated locally.
        virtual void correct();

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
