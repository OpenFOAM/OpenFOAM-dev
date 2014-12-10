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

Class
    Foam::incompressible::LESModels::DeardorffDiffStress

Group
    grpIcoLESTurbulence

Description
    Differential SGS Stress Equation Model for incompressible flows

    The DSEM uses a model version of the full balance equation for the SGS
    stress tensor to simulate the behaviour of B.
    Thus,
        \f[
            \frac{d}{dt}(B) + \nabla\cdot(U B) - \nabla\cdot(\nu_{sgs}\nabla B)
          =
            P - c_1 \frac{\epsilon}{k}B - 0.667(1 - c_1)\epsilon I
            - c_2 (P - 0.333 trP I)
        \f]

    where
    \f{eqnarray*}{
        k          &=& 0.5 tr(B)                \\
        \epsilon   &=& c_e k^{3/2}/\Delta       \\
        \epsilon/k &=& c_e k^{1/2}/\Delta       \\
        P          &=& -(B'L + L'B)             \\
        \nu_{sgs}  &=& c_k k^{1/2} \Delta       \\
        \nu_{eff}  &=& \nu_{sgs} + \nu
    \f}

SourceFiles
    DeardorffDiffStress.C

\*---------------------------------------------------------------------------*/

#ifndef DeardorffDiffStress_H
#define DeardorffDiffStress_H

#include "GenSGSStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                    Class DeardorffDiffStress Declaration
\*---------------------------------------------------------------------------*/

class DeardorffDiffStress
:
    public GenSGSStress
{
    // Private data

        dimensionedScalar ck_;
        dimensionedScalar cm_;


    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields(const volScalarField& K);

        // Disallow default bitwise copy construct and assignment
        DeardorffDiffStress(const DeardorffDiffStress&);
        DeardorffDiffStress& operator=(const DeardorffDiffStress&);


public:

    //- Runtime type information
    TypeName("DeardorffDiffStress");

    // Constructors

        //- Construct from components
        DeardorffDiffStress
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~DeardorffDiffStress()
    {}


    // Member Functions

        //- Return the effective diffusivity for B
        tmp<volScalarField> DBEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DBEff", nuSgs_ + nu())
            );
        }

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
