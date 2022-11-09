static const word divAlphaName("div(phi,alpha)");

const word alphaScheme(mesh.schemes().div(divAlphaName)[1].wordToken());

ITstream compressionScheme
(
    compressionSchemes.found(alphaScheme)
  ? mesh.schemes().div(divAlphaName)
  : ITstream
    (
        divAlphaName,
        tokenList
        {
            word("Gauss"),
            word("interfaceCompression"),
            alphaScheme,
            alphaControls.lookup<scalar>("cAlpha")
        }
    )
);
