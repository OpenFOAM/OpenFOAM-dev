/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2022 OpenFOAM Foundation
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
    Foam::functionObjects::scalarTransport

Description
    Evolves a passive scalar transport equation.

    - To specify the field name set the \c field entry
    - To employ the same numerical schemes as another field set
      the \c schemesField entry,
    - The \c diffusivity entry can be set to \c none, \c constant, \c viscosity
    - A constant diffusivity is specified with the \c D entry,
    - If a momentum transport model is available and the \c viscosity
      diffusivity option specified an effective diffusivity may be constructed
      from the laminar and turbulent viscosities using the diffusivity
      coefficients \c alphal and \c alphat:
      \verbatim
          D = alphal*nu + alphat*nut
      \endverbatim

    Example:
    \verbatim
        #includeFunc scalarTransport(T, alphal=1, alphat=1)
    \endverbatim

    For incompressible flow the passive scalar may optionally be solved with the
    MULES limiter and sub-cycling or semi-implicit in order to maintain
    boundedness, particularly if a compressive, PLIC or MPLIC convection
    scheme is used.

    Example:
    \verbatim
        #includeFunc scalarTransport(tracer, diffusion=none)

    with scheme specification:
        div(phi,tracer)     Gauss interfaceCompression vanLeer 1;

    and solver specification:
        tracer
        {
            nCorr      1;
            nSubCycles 3;

            MULESCorr       no;
            nLimiterIter    5;
            applyPrevCorr   yes;

            solver          smoothSolver;
            smoother        symGaussSeidel;
            tolerance       1e-8;
            relTol          0;

            diffusion
            {
                solver          smoothSolver;
                smoother        symGaussSeidel;
                tolerance       1e-8;
                relTol          0;
            }
        }
    \endverbatim

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    scalarTransport.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_scalarTransport_H
#define functionObjects_scalarTransport_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class scalarTransport Declaration
\*---------------------------------------------------------------------------*/

class scalarTransport
:
    public fvMeshFunctionObject
{
public:

        //- Enumeration defining the type of the diffusion
        enum class diffusionType
        {
            none,
            constant,
            viscosity
        };


private:

    // Private Data

        //- Name of field to process
        word fieldName_;

        //- Name of flux field (optional)
        word phiName_;

        //- Name of density field (optional)
        word rhoName_;

        //- Diffusion type names
        static const NamedEnum<diffusionType, 3> diffusionTypeNames_;

        //- The type of diffusion
        diffusionType diffusion_;

        //- Constant diffusion coefficient (optional)
        scalar D_;

        //- Laminar diffusion coefficient (optional)
        scalar alphal_;

        //- Turbulent diffusion coefficient (optional)
        scalar alphat_;

        //- Number of corrector iterations (optional)
        int nCorr_;

        //- Name of field whose schemes are used (optional)
        word schemesField_;

        //- The scalar field
        volScalarField s_;

        //- Switch for MULES limited solution
        bool MULES_;

        //- Stabilisation for normalisation of the interface normal
        //  needed if a compressive convection scheme is used
        const dimensionedScalar deltaN_;

        //- Scalar volumetric flux
        tmp<surfaceScalarField> tsPhi_;

        //- MULES Correction
        tmp<surfaceScalarField> tsPhiCorr0_;

        //- Switch to indicate that s has been restarted for Crank-Nicolson
        bool sRestart_;


    // Private Member Functions

        //- Return the diffusivity field
        tmp<volScalarField> D() const;

        void subCycleMULES();
        void solveMULES();


public:

    //- Runtime type information
    TypeName("scalarTransport");


    // Constructors

        //- Construct from Time and dictionary
        scalarTransport
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        scalarTransport(const scalarTransport&) = delete;


    //- Destructor
    virtual ~scalarTransport();


    // Member Functions

        //- Read the scalarTransport data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Calculate the scalarTransport
        virtual bool execute();

        //- Do nothing.
        //  The volScalarField is registered and written automatically
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const scalarTransport&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
