/*
    find the triangle in which the position is,
    when the position is known to be on the face
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool interpolationCellPointFace<Type>::findTriangle
(
    const vector& position,
    const label nFace,
    label tetPointLabels[],
    scalar phi[]
) const
{
    bool foundTriangle = false;
    vector tetPoints[3];
    const labelList& facePoints = this->pMeshFaces_[nFace];
    tetPoints[2] = this->pMeshFaceCentres_[nFace];

    label pointi = 0;

    while (pointi < facePoints.size() && !foundTriangle)
    {
        // The triangle is constructed from face center and one
        // face edge
        label nextPointLabel = (pointi + 1) % facePoints.size();

        tetPointLabels[0] = facePoints[pointi];
        tetPointLabels[1] = facePoints[nextPointLabel];

        tetPoints[0] = this->pMeshPoints_[tetPointLabels[0]];
        tetPoints[1] = this->pMeshPoints_[tetPointLabels[1]];

        vector fc = (tetPoints[0] + tetPoints[1] + tetPoints[2])/3.0;

        vector newPos = position + small*(fc-position);

        // calculate triangle edge vectors and triangle face normal
        // the 'i':th edge is opposite node i
        vector edge[3], normal[3];
        for (label i=0; i<3; i++)
        {
            label ip0 = (i+1) % 3;
            label ipp = (i+2) % 3;
            edge[i] = tetPoints[ipp]-tetPoints[ip0];
        }

        vector triangleFaceNormal = edge[1] ^ edge[2];

        // calculate edge normal (pointing inwards)
        for (label i=0; i<3; i++)
        {
            normal[i] = triangleFaceNormal ^ edge[i];
            normal[i] /= mag(normal[i]) + vSmall;
        }

        // check if position is inside triangle
        bool inside = true;
        for (label i=0; i<3; i++)
        {
            label ip = (i+1) % 3;
            inside = inside && (((newPos - tetPoints[ip]) & edge[i]) >= 0);
        }

        if (inside)
        {
            foundTriangle = true;

            // calculate phi-values
            for (label i=0; i<3; i++)
            {
                label ip = (i+1) % 3;
                scalar phiMax = max(vSmall, normal[i] & edge[ip]);
                scalar phiLength = (position-tetPoints[ip]) & normal[i];
                phi[i] = phiLength/phiMax;
            }
        }

        pointi++;
    }

    return foundTriangle;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
