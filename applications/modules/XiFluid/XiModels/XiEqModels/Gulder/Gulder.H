/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::XiEqModels::Gulder

Description
    Simple model for the equilibrium flame wrinkling \c XiEq
    based on Gulder's turbulent flame speed correlation

    References:
    \verbatim
        Gülder, Ö. L. (1991, January).
        Turbulent premixed flame propagation models
        for different combustion regimes.
        In Symposium (International) on Combustion
        (Vol. 23, No. 1, pp. 743-750). Elsevier.

        Weller, H. G. (1993).
        The development of a new flame area combustion model
        using conditional averaging.
        Thermo-fluids section report TF 9307.

        Weller, H. G., Tabor, G., Gosman, A. D., & Fureby, C. (1998, January).
        Application of a flame-wrinkling LES combustion model
        to a turbulent mixing layer.
        In Symposium (International) on combustion
        (Vol. 27, No. 1, pp. 899-907). Elsevier.
    \endverbatim

SourceFiles
    Gulder.C

\*---------------------------------------------------------------------------*/

#ifndef Gulder_H
#define Gulder_H

#include "XiEqModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{

/*---------------------------------------------------------------------------*\
                           Class Gulder Declaration
\*---------------------------------------------------------------------------*/

class Gulder
:
    public XiEqModel
{
    // Private Data

        //- Model constant
        scalar XiEqCoeff_;

        //- Schelkin effect Model constant
        scalar uPrimeCoeff_;

        //- Minimum laminar burning velocity
        const dimensionedScalar SuMin_;


protected:

        //- Update coefficients from given dictionary
        virtual bool readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("Gulder");


    // Constructors

        //- Construct from components
        Gulder
        (
            const dictionary& dict,
            const psiuMulticomponentThermo& thermo,
            const fluidThermoThermophysicalTransportModel& turbulence,
            const volScalarField& Su
        );

        //- Disallow default bitwise copy construction
        Gulder(const Gulder&) = delete;


    //- Destructor
    virtual ~Gulder();


    // Member Functions

        //- Return the flame-wrinkling XiEq
        virtual tmp<volScalarField> XiEq() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Gulder&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiEqModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
