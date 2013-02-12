#include <stdio.h>
#include <stdarg.h>
#include "grandval.h"
#include "util.h"

static int stack_size = 0;
static char **error_msg_stack = NULL;

char *make_message(const char *fmt, ...)
{
    /* Guess we need no more than 100 bytes. */
    int n, size = 100;
    char *p, *np;
    va_list ap;

    if ((p = (char *)malloc(size)) == NULL)
        return NULL;

    while (1) 
    {
        /* Try to print in the allocated space. */
        va_start(ap, fmt);
        n = vsnprintf(p, size, fmt, ap);
        va_end(ap);
        /* If that worked, return the string. */
        if (n > -1 && n < size)
            return p;
        /* Else try again with more space. */
        if (n > -1)    /* glibc 2.1 */
            size = n+1; /* precisely what is needed */
        else           /* glibc 2.0 */
            size *= 2;  /* twice the old size */
        if ((np = (char *)realloc (p, size)) == NULL) 
        {
            free(p);
            return NULL;
        } else {
            p = np;
        }
    }
}

void errmsg(const char *fmt, ...)
{
    size_t len1, len2;

    /* Guess we need no more than 100 bytes. */
    int n, size = 100;
    char *p, *np;
    va_list ap;

    if ((p = (char *)malloc(size)) == NULL)
        return;

    while (1) 
    {
        /* Try to print in the allocated space. */
        va_start(ap, fmt);
        n = vsnprintf(p, size, fmt, ap);
        va_end(ap);
        /* If that worked, return the string. */
        if (n > -1 && n < size)
            break;
        /* Else try again with more space. */
        if (n > -1)    /* glibc 2.1 */
            size = n+1; /* precisely what is needed */
        else           /* glibc 2.0 */
            size *= 2;  /* twice the old size */
        if ((np = (char *)realloc (p, size)) == NULL) 
        {
            free(p);
            return;
        } else {
            p = np;
        }
    }

    stack_size++;
    error_msg_stack = (char **)realloc(error_msg_stack, sizeof(*error_msg_stack) * stack_size);
    error_msg_stack[stack_size-1] = p;
}

void print_errmsg()
{
    int i;
    if (stack_size)
    {
        //eprintf("The following errors occured:\n");
        for (i=0; i < stack_size; i++)
        {
            eprintf("ERROR:    %s\n", error_msg_stack[i]);
            free(error_msg_stack[i]);
        }
        free(error_msg_stack);
        error_msg_stack = NULL;
        stack_size = 0;
    }
}

