#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class std::vector<std::vector<int>>+;

#pragma link C++ namespace WCP2dToy;

// #pragma link C++ class WCP2dToy::ToyEventDisplay;
// #pragma link C++ class WCP2dToy::FrameDataSource;


#pragma link C++ class WCP2dToy::ToyTiling+;
#pragma link C++ class WCP2dToy::BlobToyTiling+;
#pragma link C++ class WCP2dToy::SimpleBlobToyTiling+;
#pragma link C++ class WCP2dToy::CaveToyTiling+;
#pragma link C++ class WCP2dToy::MergeToyTiling;
#pragma link C++ class WCP2dToy::TruthToyTiling;

#pragma link C++ class WCP2dToy::ToyMatrix;
#pragma link C++ class WCP2dToy::ToyMatrixIterate;
#pragma link C++ class WCP2dToy::ToyMatrixMarkov;
#pragma link C++ class WCP2dToy::SimpleBlobToyTiling;

#endif
