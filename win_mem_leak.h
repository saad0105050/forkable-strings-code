#pragma once

#define _CRTDBG_MAP_ALLOC
#ifdef _DEBUG
#define new new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
// Replace _NORMAL_BLOCK with _CLIENT_BLOCK if you want the
// allocations to be of _CLIENT_BLOCK type
#else
#define DBG_NEW new
#endif

#define DELETE(x) if(x) delete x; x = 0;
#define DELETE_ARR(x) if(x) delete[] x; x = 0;

#include <stdlib.h>  
#include <crtdbg.h>  
