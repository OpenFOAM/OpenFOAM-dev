/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

    Foam::
    compressible::
    filmPyrolysisRadiativeCoupledMixedFvPatchScalarField

Description
    Mixed boundary condition for temperature, to be used in the flow and
    pyrolysis regions when a film region model is used.


    Example usage:
        myInterfacePatchName
        {
            type            filmPyrolysisRadiativeCoupledMixed;
            Tnbr            T;
            kappa           fluidThermo;
            Qr              Qr;
            kappaName       none;
            filmDeltaDry    0.0;
            filmDeltaWet    3e-4;
            value           $internalField;
        }

    Needs to be on underlying mapped(Wall)FvPatch.
    It calculates local field as
         ratio = (filmDelta - filmDeltaDry)/(filmDeltaWet - filmDeltaDry),

    when ratio = 1 is considered wet and the film temperarture is fixed at
    the wall. If ratio = 0 (dry) it emulates the normal radiative solid BC.

    In between ratio 0 and 1 the gradient and value contributions are
    weighted using the ratio field in the followig way:

    qConv = ratio*htcwfilm*(Tfilm - *this);
    qRad = (1.0 - ratio)*Qr;

    Then the solid can gain or loose energy through radiation or conduction
    towards the film.

    Note: kappa : heat conduction at patch.
        Gets supplied how to lookup/calculate kappa:
    - 'lookup' : lookup volScalarField (or volSymmTensorField) with name
    - 'basicThermo' : use basicThermo and compressible::RASmodel to calculate K
    - 'solidThermo' : use basicSolidThermo K()

    Qr is the radiative flux defined in the radiation model.


SourceFiles
    filmPyrolysisRadiativeCoupledMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef filmPyrolysisRadiativeCoupledMixedFvPatchScalarField_H
#define filmPyrolysisRadiativeCoupledMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "thermoSingleLayer.H"
#include "pyrolysisModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class filmPyrolysisRadiativeCoupledMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase
{
public:

    typedef Foam::regionModels::surfaceFilmModels::thermoSingleLayer
        filmModelType;

    typedef Foam::regionModels::pyrolysisModels::pyrolysisModel
        pyrolysisModelType;


private:

    // Private data

        //- Name of film region
        const word filmRegionName_;

        //- Name of pyrolysis region
        const word pyrolysisRegionName_;

        //- Name of field on the neighbour region
        const word TnbrName_;

        //- Name of the radiative heat flux
        const word QrName_;

        //- Convective Scaling Factor (as determined by Prateep's tests)
        const scalar convectiveScaling_;

        //- Minimum delta film to be consired dry
        const scalar filmDeltaDry_;

        //- Maximum delta film to be consired wet
        const scalar filmDeltaWet_;

        //- Retrieve film model from the database
        const filmModelType& filmModel() const;

        //- Retrieve pyrolysis model from the database
        const pyrolysisModelType& pyrModel() const;


public:

    //- Runtime type information
    TypeName("filmPyrolysisRadiativeCoupledMixed");


    // Constructors

        //- Construct from patch and internal field
        filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
        (
            const
            filmPyrolysisRadiativeCoupledMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
        (
            const filmPyrolysisRadiativeCoupledMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Get corresponding K field
        tmp<scalarField> K() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
