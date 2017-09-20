/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::IntegrationScheme

Description
    Base for a set of schemes which integrate simple ODEs which arise from
    semi-implcit rate expressions.

        \f[
            \frac{d \phi}{d t} = A - B \phi
        \f]

    The methods are defined in terms of the factor \f$f\f$ by which the result
    of a trivial explicit integration is multipled.

        \f[
            \Delta \phi = (A - B \phi^n) f(\Delta t, B)
        \f]

    This class also facilitates integration in stages. If the explicit and
    implicit coefficients, \f$A\f$ and \f$B\f$, are a summation of differing
    contributions, \f$\sum \alpha_i\f$ and \f$\sum \beta_i\f$, then the
    integration can be split up to detemine the effect of each contribution.

        \f[
            \frac{d \phi}{d t} = \left( \sum \alpha_i \right) -
            \left( \sum \beta_i \right) \phi
        \f]
        \f[
            \Delta \phi_i = (\alpha_i - \beta_i \phi^n) f(\Delta t, B)
        \f]

    Note that the sum of implicit coefficients \f$B\f$ is needed for each stage
    of the integral, whilst the sum of explicit coefficients \f$A\f$ is not.

SourceFiles
    IntegrationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef IntegrationScheme_H
#define IntegrationScheme_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class IntegrationScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class IntegrationScheme
{
public:

    //- Runtime type information
    TypeName("integrationScheme");


    //- Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             IntegrationScheme,
             word,
             (),
             ()
         );


    // Constructors

        //- Construct
        IntegrationScheme();

        //- Construct and return clone
        virtual autoPtr<IntegrationScheme<Type>> clone() const = 0;


    // Selectors

        //- Select an integration scheme
        static autoPtr<IntegrationScheme> New
        (
            const word& phiName,
            const dictionary& dict
        );


    //- Destructor
    virtual ~IntegrationScheme();


    // Member Functions

        //- Perform the full integration
        Type delta
        (
            const Type& phi,
            const scalar dt,
            const Type& Alpha,
            const scalar Beta
        ) const;

        //- Perform a stage of the integration
        Type delta
        (
            const Type& phi,
            const scalar dt,
            const scalar Beta,
            const Type& alpha,
            const scalar beta
        ) const;

        //- Return the factor
        virtual scalar factor(const scalar dt, const scalar Beta) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeIntegrationScheme(Type)                                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(IntegrationScheme<Type>, 0);           \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        IntegrationScheme<Type>,                                               \
        word                                                                   \
    );


#define makeIntegrationSchemeType(SS, Type)                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                          \
                                                                               \
    IntegrationScheme<Type>::addwordConstructorToTable<SS<Type>>               \
        add##SS##Type##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IntegrationScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
