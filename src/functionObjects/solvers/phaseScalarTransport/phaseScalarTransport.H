/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2021 OpenFOAM Foundation
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
    Foam::functionObjects::phaseScalarTransport

Description
    Evolves a passive scalar transport equation within one phase of a
    multiphase simulation. The scalar is considered to be a phase-intensive
    property; i.e., its value represents an amount per-unit of the phase. In
    addition to the scalar, the function also writes out the product of the
    volume fraction and the scalar, as this provides a phase-extensive field
    which is often more convenient to post-process.

    Most entries are the same as for the \c scalarTransport function. Refer to
    its documentation for details. Entries specific to this function are
    detailed below. Note that the phase-name will be determined by stripping
    the extension from the supplied field name.

    If \c alphaPhi is specified and found in the database then it will be used
    to transport the field. This is the likely mode of operation if this
    function is used with Euler-Euler solvers, in which \c alphaPhi -like
    fluxes are available. If \c alphaPhi is not found, then a pressure-like
    equation will be solved in order to construct it so that it exactly matches
    the time-derivative of the phase volume or mass. This is likely to be
    necessary in volume-of-fluid solvers where \c alphaPhi is not part of the
    solution procedure. The pressure field name will be required in this case.

Usage
    \table
        Property        | Description                        | Req'd? | Default
        alpha           | Name of the volume-fraction field  | no\\
                        | alpha.<phase-name>
        alphaPhi        | Name of the phase-flux field       | no\\
                        | alphaPhi.<phase-name>
        p               | Name of the pressure field         | no     | p
        residualAlpha   | Small volume fraction used to stabilise the solution\\
                        | no | rootSmall
        writeAlphaField | Also write out alpha multiplied by the field\\
                        | no | true
    \endtable

    Example specification for interFoam:
    \verbatim
    phaseScalarTransport1
    {
        type            phaseScalarTransport;
        libs            ("libsolverFunctionObjects.so");

        field           s.water;
        p               p_rgh;
    }
    \endverbatim

    Example specification for multiphaseEulerFoam:
    \verbatim
    phaseScalarTransport1
    {
        type            phaseScalarTransport;
        libs            ("libsolverFunctionObjects.so");

        field           s.water;
        alphaPhi        alphaRhoPhi.water;
        rho             thermo:rho.water;
    }
    \endverbatim

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::scalarTransport

SourceFiles
    phaseScalarTransport.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_phaseScalarTransport_H
#define functionObjects_phaseScalarTransport_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class phaseScalarTransport Declaration
\*---------------------------------------------------------------------------*/

class phaseScalarTransport
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Name of field to process
        const word fieldName_;

        //- Name of the phase in which to solve
        const word phaseName_;

        //- Name of phase volume-fraction field (optional)
        word alphaName_;

        //- Name of the phase flux field (optional)
        word alphaPhiName_;

        //- Name of the mixture flux field (optional)
        word phiName_;

        //- Name of density field (optional)
        word rhoName_;

        //- Name of the pressure field (optional)
        word pName_;

        //- Diffusion coefficient (optional)
        scalar D_;

        //- Flag to indicate whether a constant, uniform D_ is specified
        bool constantD_;

        //- Laminar diffusion coefficient (optional)
        scalar alphaD_;

        //- Turbulent diffusion coefficient (optional)
        scalar alphaDt_;

        //- Number of corrector iterations (optional)
        int nCorr_;

        //- Residual volume-fraction
        scalar residualAlpha_;

        //- Name of field whose schemes are used (optional)
        word schemesField_;

        //- Flag to indicate whether to write the field multiplied by the phase
        //  fraction
        bool writeAlphaField_;

        //- The field
        volScalarField s_;

        //- The field multiplied by the phase fraction
        autoPtr<volScalarField> alphaSPtr_;

        //- Potential field used to generate the phase flux
        autoPtr<volScalarField> PhiPtr_;


    // Private Member Functions

        //- Return the potential field used to generate the phase flux.
        //  Constructed on demand.
        volScalarField& Phi();

        //- Return the phase flux. Tries to look it up, and generates it if the
        //  lookup fails. The generation creates a crude guess for alphaPhi
        //  then solves a Laplacian to correct the flux to match the time
        //  derivative of alpha.
        tmp<surfaceScalarField> alphaPhi();

        //- Return the diffusivity field
        tmp<volScalarField> D(const surfaceScalarField& alphaPhi) const;


public:

    //- Runtime type information
    TypeName("phaseScalarTransport");


    // Constructors

        //- Construct from Time and dictionary
        phaseScalarTransport
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~phaseScalarTransport();


    // Member Functions

        //- Read the settings from the given dictionary
        virtual bool read(const dictionary&);

        //- Solve for the evolution of the field
        virtual bool execute();

        //- Do nothing. The field is registered and written automatically.
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
