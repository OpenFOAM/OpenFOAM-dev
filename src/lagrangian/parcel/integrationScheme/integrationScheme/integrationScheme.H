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
    Foam::integrationScheme

Description
    Base for a set of schemes which integrate simple ODEs which arise from
    semi-implicit rate expressions.

        \f[
            \frac{d \phi}{d t} = A - B \phi
        \f]

    The methods are defined in terms of the effective time-step \f$\Delta
    t_e\f$ by which the explicit rate is multiplied. The effective time-step is
    a function of the actual time step and the implicit coefficient, which must
    be implemented in each derived scheme.

        \f[
            \Delta t_e = f(\Delta t, B)
        \f]
        \f[
            \Delta \phi = (A - B \phi^n) \Delta t_e
        \f]

    This class also facilitates integration in stages. If the explicit and
    implicit coefficients, \f$A\f$ and \f$B\f$, are a summation of differing
    contributions, \f$\sum \alpha_i\f$ and \f$\sum \beta_i\f$, then the
    integration can be split up to determine the effect of each contribution.

        \f[
            \frac{d \phi_i}{d t} = \alpha_i - \beta_i \phi
        \f]
        \f[
            \Delta \phi_i = \alpha_i \Delta t -
            \beta_i \int_0^{\Delta t} \phi d t
        \f]
        \f[
            \Delta \phi_i = (\alpha_i - \beta_i \phi^n) \Delta t -
            \beta_i (A - B \phi^n) \int_0^{\Delta t} t_e dt
        \f]

    These partial calculations are defined in terms of the integral of the
    effective time-step, \f$\int_0^{\Delta t} t_e dt\f$, which is also
    implemented in every derivation.

SourceFiles
    integrationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef integrationScheme_H
#define integrationScheme_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class integrationScheme Declaration
\*---------------------------------------------------------------------------*/

class integrationScheme
{
public:

    //- Runtime type information
    TypeName("integrationScheme");


    //- Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             integrationScheme,
             word,
             (),
             ()
         );


    // Constructors

        //- Construct
        integrationScheme();

        //- Construct and return clone
        virtual autoPtr<integrationScheme> clone() const = 0;


    // Selectors

        //- Select an integration scheme
        static autoPtr<integrationScheme> New
        (
            const word& phiName,
            const dictionary& dict
        );


    //- Destructor
    virtual ~integrationScheme();


    // Member Functions

        //- Perform the integration explicitly
        template<class Type>
        static Type explicitDelta
        (
            const Type& phi,
            const scalar dtEff,
            const Type& Alpha,
            const scalar Beta
        );

        //- Perform the integration
        template<class Type>
        Type delta
        (
            const Type& phi,
            const scalar dt,
            const Type& Alpha,
            const scalar Beta
        ) const;

        //- Perform a part of the integration
        template<class Type>
        Type partialDelta
        (
            const Type& phi,
            const scalar dt,
            const Type& Alpha,
            const scalar Beta,
            const Type& alphai,
            const scalar betai
        ) const;

        //- Return the integration effective time step
        virtual scalar dtEff(const scalar dt, const scalar Beta) const = 0;

        //- Return the integral of the effective time step
        virtual scalar sumDtEff(const scalar dt, const scalar Beta) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "integrationSchemeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
