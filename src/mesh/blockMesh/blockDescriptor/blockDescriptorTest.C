namespace Foam
{
    inline point faceCorr(const point& p)
    {
        return p/mag(p);
    }
}


Foam::label Foam::blockDescriptor::correctFacePoints
(
    FixedList<List<point>, 6>& facePoints
) const
{
    forAll(facePoints, facei)
    {
        forAll(facePoints[facei], pointi)
        {
            facePoints[facei][pointi] = faceCorr(facePoints[facei][pointi]);
        }
    }

    return 6;
}
