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
    Foam::XiModels::transport

Description
    Transport model for flame wrinkling \c Xi

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
    transport.C

\*---------------------------------------------------------------------------*/

#ifndef transport_H
#define transport_H

#include "XiModel.H"
#include "XiEqModel.H"
#include "XiGModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiModels
{

/*---------------------------------------------------------------------------*\
                          Class transport Declaration
\*---------------------------------------------------------------------------*/

class transport
:
    public XiModel
{
    // Private Data

        scalar XiShapeCoeff_;

        autoPtr<XiEqModel> XiEqModel_;
        autoPtr<XiGModel> XiGModel_;


protected:

        //- Update coefficients from given dictionary
        virtual bool readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("transport");


    // Constructors

        //- Construct from components
        transport
        (
            const dictionary& dict,
            const psiuMulticomponentThermo& thermo,
            const fluidThermoThermophysicalTransportModel& turbulence,
            const volScalarField& Su
        );

        //- Disallow default bitwise copy construction
        transport(const transport&) = delete;


    //- Destructor
    virtual ~transport();


    // Member Functions

        //- Return the flame diffusivity
        virtual tmp<volScalarField> Db() const;

        //- Add Xi to the multivariateSurfaceInterpolationScheme table
        virtual void addXi
        (
            multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
        )
        {
            fields.add(Xi_);
        }

        //- Correct the flame-wrinkling Xi
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const transport&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
