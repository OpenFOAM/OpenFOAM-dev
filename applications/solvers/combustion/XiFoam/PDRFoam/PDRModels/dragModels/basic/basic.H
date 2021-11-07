/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::PDRDragModels::basic

Description
    Basic sub-grid obstacle drag model.
    Details supplied by J Puttock 2/7/06.

    <b> Sub-grid drag term </b>

    The resistance term (force per unit of volume) is given by:

    \f[
        R = -\frac{1}{2} \rho \vert \dwea{\vec{U}} \vert \dwea{\vec{U}}.D
    \f]

    where:

        \f$ D \f$ is the tensor field "CR" in \f$ m^{-1} \f$

        This is term is treated implicitly in UEqn.H

    <b> Sub-grid turbulence generation </b>

    The turbulence source term \f$ G_{R} \f$ occurring in the
    \f$ \kappa-\epsilon \f$ equations for the generation of turbulence due
    to interaction with unresolved obstacles :

    \f$ G_{R} = C_{s}\beta_{\nu}
    \mu_{eff} A_{w}^{2}(\dwea{\vec{U}}-\dwea{\vec{U}_{s}})^2 + \frac{1}{2}
    \rho \vert \dwea{\vec{U}} \vert \dwea{\vec{U}}.T.\dwea{\vec{U}} \f$

    where:

        \f$ C_{s} \f$ = 1

        \f$ \beta_{\nu} \f$ is the volume porosity (file "betav").

        \f$ \mu_{eff} \f$ is the effective viscosity.

        \f$ A_{w}^{2}\f$ is the obstacle surface area per unit of volume
        (file "Aw").

        \f$ \dwea{\vec{U}_{s}} \f$ is the slip velocity and is considered
        \f$ \frac{1}{2}. \dwea{\vec{U}} \f$.

        \f$ T \f$ is a tensor in the file CT.

    The term \f$ G_{R} \f$ is treated explicitly in the \f$ \kappa-\epsilon
    \f$ Eqs in the \link PDRkEpsilon.C \endlink file.


SourceFiles
    basic.C

\*---------------------------------------------------------------------------*/

#ifndef basic_H
#define basic_H

#include "PDRDragModel.H"
#include "XiEqModel.H"


namespace Foam
{
namespace PDRDragModels
{

/*---------------------------------------------------------------------------*\
                            Class basic Declaration
\*---------------------------------------------------------------------------*/

class basic
:
    public PDRDragModel
{
    // Private Data

        dimensionedScalar Csu;
        dimensionedScalar Csk;

        volScalarField Aw_;
        volSymmTensorField CR_;


public:

    //- Runtime type information
    TypeName("basic");


    // Constructors

        //- Construct from components
        basic
        (
            const dictionary& PDRProperties,
            const compressible::RASModel& turbulence,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi
        );

        //- Disallow default bitwise copy construction
        basic(const basic&);


    //- Destructor
    virtual ~basic();


    // Member Functions

        //- Return the momentum drag coefficient
        virtual tmp<volSymmTensorField> Dcu() const;

        //- Return the momentum drag turbulence generation rate
        virtual tmp<volScalarField> Gk() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& PDRProperties);

        //- Write fields
        void writeFields() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const basic&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PDRDragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
