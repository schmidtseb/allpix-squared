#include "vector"
#ifdef __CINT__ 
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ class vector<allpix::PixelHit*>+;
#pragma link C++ class vector<allpix::PixelHit*>::*;
#ifdef G__VECTOR_HAS_CLASS_ITERATOR
#pragma link C++ operators vector<allpix::PixelHit*>::iterator;
#pragma link C++ operators vector<allpix::PixelHit*>::const_iterator;
#pragma link C++ operators vector<allpix::PixelHit*>::reverse_iterator;
#endif
#endif
