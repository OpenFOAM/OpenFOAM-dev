/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "MixedDiffuseSpecular.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MixedDiffuseSpecular<CloudType>::MixedDiffuseSpecular
(
    const dictionary& dict,
    CloudType& cloud
)
:
    WallInteractionModel<CloudType>(dict, cloud, typeName),
    diffuseFraction_(readScalar(this->coeffDict().lookup("diffuseFraction")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MixedDiffuseSpecular<CloudType>::~MixedDiffuseSpecular()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::MixedDiffuseSpecular<CloudType>::correct
(
    typename CloudType::parcelType& p
)
{
    vector& U = p.U();

    scalar& Ei = p.Ei();

    label typeId = p.typeId();

    const label wppIndex = p.patch();

    const polyPatch& wpp = p.mesh().boundaryMesh()[wppIndex];

    label wppLocalFace = wpp.whichFace(p.face());

    const vector nw = p.normal();

    // Normal velocity magnitude
    scalar U_dot_nw = U & nw;

    CloudType& cloud(this->owner());

    Random& rndGen(cloud.rndGen());

    if (diffuseFraction_ > rndGen.scalar01())
    {
        // Diffuse reflection

        // Wall tangential velocity (flow direction)
        vector Ut = U - U_dot_nw*nw;

        while (mag(Ut) < small)
        {
            // If the incident velocity is parallel to the face normal, no
            // tangential direction can be chosen.  Add a perturbation to the
            // incoming velocity and recalculate.

            U = vector
            (
                U.x()*(0.8 + 0.2*rndGen.scalar01()),
                U.y()*(0.8 + 0.2*rndGen.scalar01()),
                U.z()*(0.8 + 0.2*rndGen.scalar01())
            );

            U_dot_nw = U & nw;

            Ut = U - U_dot_nw*nw;
        }

        // Wall tangential unit vector
        vector tw1 = Ut/mag(Ut);

        // Other tangential unit vector
        vector tw2 = nw^tw1;

        scalar T = cloud.boundaryT().boundaryField()[wppIndex][wppLocalFace];

        scalar mass = cloud.constProps(typeId).mass();

        direction iDof = cloud.constProps(typeId).internalDegreesOfFreedom();

        U =
            sqrt(physicoChemical::k.value()*T/mass)
           *(
                rndGen.GaussNormal()*tw1
              + rndGen.GaussNormal()*tw2
              - sqrt(-2.0*log(max(1 - rndGen.scalar01(), vSmall)))*nw
            );

        U += cloud.boundaryU().boundaryField()[wppIndex][wppLocalFace];

        Ei = cloud.equipartitionInternalEnergy(T, iDof);
    }
    else
    {
        // Specular reflection

        if (U_dot_nw > 0.0)
        {
            U -= 2.0*U_dot_nw*nw;
        }
    }

}


// ************************************************************************* //
