transformer transforms;

{
    wordReList simpleTransformations;
    List<Tuple2<word, string>> transformations;
    dictArgList(transformationString, simpleTransformations, transformations);

    forAll(transformations, i)
    {
        if (transformations[i].first() == "translate")
        {
            const vector v(IStringStream(transformations[i].second())());
            transforms = transformer::translation(v) & transforms;
        }
        else if (transformations[i].first() == "rotate")
        {
            Pair<vector> n1n2(IStringStream(transformations[i].second())());

            n1n2[0] /= mag(n1n2[0]);
            n1n2[1] /= mag(n1n2[1]);

            transforms =
                transformer::rotation(rotationTensor(n1n2[0], n1n2[1]))
              & transforms;
        }
        else if (transformations[i].first() == "Rx")
        {
            const scalar a
            (
                readScalar(IStringStream(transformations[i].second())())
            );
            transforms = transformer::rotation(Rx(degToRad(a))) & transforms;
        }
        else if (transformations[i].first() == "Ry")
        {
            const scalar a
            (
                readScalar(IStringStream(transformations[i].second())())
            );
            transforms = transformer::rotation(Ry(degToRad(a))) & transforms;
        }
        else if (transformations[i].first() == "Rz")
        {
            const scalar a
            (
                readScalar(IStringStream(transformations[i].second())())
            );
            transforms = transformer::rotation(Rz(degToRad(a))) & transforms;
        }
        else if (transformations[i].first() == "Ra")
        {
            IStringStream istr(transformations[i].second());
            const vector v(istr);
            const scalar a(readScalar(istr));
            transforms = transformer::rotation(Ra(v, degToRad(a))) & transforms;
        }
        else if (transformations[i].first() == "scale")
        {
            const vector v(IStringStream(transformations[i].second())());
            transforms =
                transformer::scaling(diagTensor(v.x(), v.y(), v.z()))
              & transforms;
        }
        else
        {
            args.printUsage();
            FatalErrorInFunction
                << "Unknown transformation " << transformations[i].first()
                << exit(FatalError);
        }
    }
}
