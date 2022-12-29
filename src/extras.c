/**
 * @file extras.c
 *
 * See https://github.com/emscripten-core/emscripten/blob/main/test/mallinfo.cpp
 */

#include <unistd.h>
#include <emscripten.h>

typedef struct
{
	int arena;    /* non-mmapped space allocated from system */
	int ordblks;  /* number of free chunks */
	int smblks;   /* always 0 */
	int hblks;    /* always 0 */
	int hblkhd;   /* space in mmapped regions */
	int usmblks;  /* maximum total allocated space */
	int fsmblks;  /* always 0 */
	int uordblks; /* total allocated space */
	int fordblks; /* total free space */
	int keepcost; /* releasable (via malloc_trim) space */
}
s_mallinfo;

extern s_mallinfo mallinfo();

EMSCRIPTEN_KEEPALIVE
unsigned int getUsedMemory()
{
    s_mallinfo info = mallinfo();
    unsigned int dynamicTop = (unsigned int) sbrk(0);
    return dynamicTop - info.fordblks;
}

/* vim: set ft=c.doxygen: */
