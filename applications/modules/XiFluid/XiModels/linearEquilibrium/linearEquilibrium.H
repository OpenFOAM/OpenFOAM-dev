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
    Foam::XiModels::linearEquilibrium

Description
    Simple equilibrium model for flame wrinkling \c Xi
    with a linear profile function to give a plausible distribution.

    Note that this model does not guarantee that the effective turbulent flame
    speed will correspond exactly to the equilibrium \c Xi because \c Xi will
    only equal the equilibrium value where \c b = 0.5.  It would be simple to
    implement alternative profiles where \c Xi equals the equilibrium value over
    a wider range of \c b within the flame-brush which should ensure the
    effective flame-speed corresponds to the equilibrium \c Xi more closely.

    See Foam::XiModel for more details on flame wrinkling modelling.

    References:
    \verbatim
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
    linearEquilibrium.C

\*---------------------------------------------------------------------------*/

#ifndef linearEquilibrium_H
#define linearEquilibrium_H

#include "equilibrium.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiModels
{

/*---------------------------------------------------------------------------*\
                          Class linearEquilibrium Declaration
\*---------------------------------------------------------------------------*/

class linearEquilibrium
:
    public equilibrium
{
    // Private Data

        scalar XiShapeCoeff_;


protected:

        //- Update coefficients from given dictionary
        virtual bool readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("linearEquilibrium");


    // Constructors

        //- Construct from components
        linearEquilibrium
        (
            const dictionary& dict,
            const psiuMulticomponentThermo& thermo,
            const fluidThermoThermophysicalTransportModel& turbulence,
            const volScalarField& Su
        );

        //- Disallow default bitwise copy construction
        linearEquilibrium(const linearEquilibrium&) = delete;


    //- Destructor
    virtual ~linearEquilibrium();


    // Member Functions

        //- Correct the flame-wrinkling Xi
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const linearEquilibrium&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
