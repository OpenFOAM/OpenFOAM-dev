/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2024 OpenFOAM Foundation
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
    Foam::momentumTransportModel

Description
    Abstract base class for turbulence models (RAS, LES and laminar).

SourceFiles
    momentumTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef momentumTransportModel_H
#define momentumTransportModel_H

#include "IOdictionary.H"
#include "primitiveFieldsFwd.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "geometricOneField.H"
#include "viscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                   Class momentumTransportModel Declaration
\*---------------------------------------------------------------------------*/

class momentumTransportModel
:
    public IOdictionary
{
protected:

    // Protected data

        const Time& runTime_;
        const fvMesh& mesh_;

        const volVectorField& U_;
        const surfaceScalarField& alphaRhoPhi_;
        const surfaceScalarField& phi_;

        const viscosity& viscosity_;


    // Protected member functions

        static typeIOobject<IOdictionary> readModelDict
        (
            const objectRegistry& obr,
            const word& group,
            bool registerObject = false
        );


        template<class MomentumTransportModel>
        static inline autoPtr<MomentumTransportModel> New
        (
            const typename MomentumTransportModel::alphaField& alpha,
            const typename MomentumTransportModel::rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );


public:

    //- Runtime type information
    TypeName("momentumTransport");


    // Constructors

        //- Construct from components
        momentumTransportModel
        (
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );

        //- Disallow default bitwise copy construction
        momentumTransportModel(const momentumTransportModel&) = delete;


    //- Destructor
    virtual ~momentumTransportModel()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read() = 0;

        const Time& time() const
        {
            return runTime_;
        }

        const fvMesh& mesh() const
        {
            return mesh_;
        }

        inline word groupName(const word& name) const
        {
            return IOobject::groupName(name, alphaRhoPhi_.group());
        }

        //- Helper function to return the name of the turbulence G field
        inline word GName() const
        {
            return typedName("G");
        }

        //- Access function to velocity field
        inline const volVectorField& U() const
        {
            return U_;
        }

        //- Access function to phase flux field
        inline const surfaceScalarField& alphaRhoPhi() const
        {
            return alphaRhoPhi_;
        }

        //- Return the volumetric flux field
        virtual tmp<surfaceScalarField> phi() const;

        //- Access function to fluid properties
        const class Foam::viscosity& properties() const
        {
            return viscosity_;
        }

        //- Return the wall distance
        const volScalarField& y() const;

        //- Return the near wall distance
        const volScalarField::Boundary& yb() const;

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return this->viscosity_.nu();
        }

        //- Return the laminar viscosity on patchi
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return this->viscosity_.nu(patchi);
        }

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const = 0;

        //- Return the turbulence viscosity on patch
        virtual tmp<scalarField> nut(const label patchi) const = 0;

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const = 0;

        //- Return the effective viscosity on patch
        virtual tmp<scalarField> nuEff(const label patchi) const = 0;

        //- Return the Reynolds stress tensor [m^2/s^2]
        virtual tmp<volSymmTensorField> R() const = 0;

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const = 0;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const = 0;

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const = 0;

        //- Return the Reynolds stress tensor [m^2/s^2]

        //- Return the explicit stress correction matrix
        //  for the momentum equation.
        //  Combined with the implicit part by divDevTau
        //  Also creates and registers the devTauCorrFlux field
        //  which is cached in the matrix if fluxRequired.
        virtual tmp<fvVectorMatrix> divDevTauCorr
        (
            const tmp<surfaceVectorField>& devTauCorr,
            volVectorField& U
        ) const;

        //- Validate the fields after construction
        //  Update derived fields as required
        virtual void validate();

        //- Predict the momentum transport coefficients if possible
        //  without solving momentum transport model equations
        virtual void predict() = 0;

        //- Solve the momentum transport model equations
        //  and correct the momentum transport coefficients
        virtual void correct() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const momentumTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "momentumTransportModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
