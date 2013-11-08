//

#include <stdlib.h>

// vector of 32-bit intergers (added for 64-bit portability)
struct veci_t 
{
    int    size;
    int    cap;
    int*   ptr;
};

typedef struct veci_t veci;

static inline void veci_new (veci* v) 
{
    v->size = 0;
    v->cap  = 4;
    v->ptr  = (int*)malloc(sizeof(int)*v->cap);
}

static inline void   veci_delete (veci* v)          { free(v->ptr);   }
static inline int*   veci_begin  (veci* v)          { return v->ptr;  }
static inline int    veci_size   (veci* v)          { return v->size; }
static inline void   veci_resize (veci* v, int k)   { v->size = k;    } // only safe to shrink !!
static inline void   veci_push   (veci* v, int e)
{
    if (v->size == v->cap) {
        int newsize = v->cap * 2+1;
        v->ptr = (int*)realloc(v->ptr,sizeof(int)*newsize);
        v->cap = newsize; }
    v->ptr[v->size++] = e;
}


// vector of 32- or 64-bit pointers
struct vecp_t 
{
    int    size;
    int    cap;
    void** ptr;
};

typedef struct vecp_t vecp;

static inline void vecp_new (vecp* v) {
    v->size = 0;
    v->cap  = 4;
    v->ptr  = (void**)malloc(sizeof(void*)*v->cap);
}

static inline void   vecp_delete (vecp* v)          { free(v->ptr);   }
static inline void** vecp_begin  (vecp* v)          { return v->ptr;  }
static inline int    vecp_size   (vecp* v)          { return v->size; }
static inline void   vecp_resize (vecp* v, int   k) { v->size = k;    } // only safe to shrink !!
static inline void   vecp_push   (vecp* v, void* e)
{
    if (v->size == v->cap) {
        int newsize = v->cap * 2+1;
        v->ptr = (void**)realloc(v->ptr,sizeof(void*)*newsize);
        v->cap = newsize; }
    v->ptr[v->size++] = e;
}
