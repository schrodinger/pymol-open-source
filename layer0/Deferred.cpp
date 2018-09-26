#include "Deferred.h"

#include "MemoryDebug.h"

void CDeferred::exec(){
    if(fn){
        fn(this);
    }
}
