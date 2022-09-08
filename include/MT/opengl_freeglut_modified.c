/*
 * freeglut_main.c
 *
 * The windows message processing methods.
 *
 * Copyright (c) 1999-2000 Pawel W. Olszta. All Rights Reserved.
 * Written by Pawel W. Olszta, <olszta@sourceforge.net>
 * Creation date: Fri Dec 3 1999
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * PAWEL W. OLSZTA BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#define XLIB_ILLEGAL_ACCESS
#include <GL/freeglut.h>
#include "freeglut_internal.h"
#include <errno.h>
#include <stdarg.h>
#if TARGET_HOST_WIN32
#    define VFPRINTF(s,f,a) vfprintf((s),(f),(a))
#else
#    if HAVE_VPRINTF
#        define VFPRINTF(s,f,a) vfprintf((s),(f),(a))
#    elif HAVE_DOPRNT
#        define VFPRINTF(s,f,a) _doprnt((f),(a),(s))
#    else
#        define VFPRINTF(s,f,a)
#    endif
#endif

#if TARGET_HOST_WINCE

typedef struct GXDisplayProperties GXDisplayProperties;
typedef struct GXKeyList GXKeyList;
#include <gx.h>

typedef struct GXKeyList (*GXGETDEFAULTKEYS)(int);
typedef int (*GXOPENINPUT)();

GXGETDEFAULTKEYS GXGetDefaultKeys_ = NULL;
GXOPENINPUT GXOpenInput_ = NULL;

struct GXKeyList gxKeyList;

#endif

/*
 * Try to get the maximum value allowed for ints, falling back to the minimum
 * guaranteed by ISO C99 if there is no suitable header.
 */
#if HAVE_LIMITS_H
#    include <limits.h>
#endif
#ifndef INT_MAX
#    define INT_MAX 32767
#endif

#ifndef MIN
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#endif




/*
 * Indicates whether Joystick events are being used by ANY window.
 *
 * The current mechanism is to walk all of the windows and ask if
 * there is a joystick callback.  We have a short-circuit early
 * return if we find any joystick handler registered.
 *
 * The real way to do this is to make use of the glutTimer() API
 * to more cleanly re-implement the joystick API.  Then, this code
 * and all other "joystick timer" code can be yanked.
 *
 */
static void fghCheckJoystickCallback( SFG_Window* w, SFG_Enumerator* e)
{
    if( FETCH_WCB( *w, Joystick ) )
    {
        e->found = GL_TRUE;
        e->data = w;
    }
    fgEnumSubWindows( w, fghCheckJoystickCallback, e );
}
static int fghHaveJoystick( void )
{
    SFG_Enumerator enumerator;

    enumerator.found = GL_FALSE;
    enumerator.data = NULL;
    fgEnumWindows( fghCheckJoystickCallback, &enumerator );
    return !!enumerator.data;
}
static void fghHavePendingRedisplaysCallback( SFG_Window* w, SFG_Enumerator* e)
{
    if( w->State.Redisplay )
    {
        e->found = GL_TRUE;
        e->data = w;
    }
    fgEnumSubWindows( w, fghHavePendingRedisplaysCallback, e );
}
static int fghHavePendingRedisplays (void)
{
    SFG_Enumerator enumerator;

    enumerator.found = GL_FALSE;
    enumerator.data = NULL;
    fgEnumWindows( fghHavePendingRedisplaysCallback, &enumerator );
    return !!enumerator.data;
}
/*
 * Returns the number of GLUT ticks (milliseconds) till the next timer event.
 */
static long fghNextTimer( void )
{
    long ret = INT_MAX;
    SFG_Timer *timer = (struct tagSFG_Timer*)fgState.Timers.First;

    if( timer )
        ret = timer->TriggerTime - fgElapsedTime();
    if( ret < 0 )
        ret = 0;

    return ret;
}
/*
 * Does the magic required to relinquish the CPU until something interesting
 * happens.
 */
static void fghSleepForEvents( void )
{
    long msec;

    if( fgState.IdleCallback || fghHavePendingRedisplays( ) )
        return;

    msec = fghNextTimer( );
    /* XXX Use GLUT timers for joysticks... */
    /* XXX Dumb; forces granularity to .01sec */
    if( fghHaveJoystick( ) && ( msec > 10 ) )     
        msec = 10;

#if TARGET_HOST_UNIX_X11
    /*
     * Possibly due to aggressive use of XFlush() and friends,
     * it is possible to have our socket drained but still have
     * unprocessed events.  (Or, this may just be normal with
     * X, anyway?)  We do non-trivial processing of X events
     * after the event-reading loop, in any case, so we
     * need to allow that we may have an empty socket but non-
     * empty event queue.
     */
    if( ! XPending( fgDisplay.Display ) )
    {
        fd_set fdset;
        int err;
        int socket;
        struct timeval wait;

        socket = ConnectionNumber( fgDisplay.Display );
        FD_ZERO( &fdset );
        FD_SET( socket, &fdset );
        wait.tv_sec = msec / 1000;
        wait.tv_usec = (msec % 1000) * 1000;
        err = select( socket+1, &fdset, NULL, NULL, &wait );

        if( ( -1 == err ) && ( errno != EINTR ) )
            fgWarning ( "freeglut select() error: %d", errno );
    }
#elif TARGET_HOST_WIN32 || TARGET_HOST_WINCE
    MsgWaitForMultipleObjects( 0, NULL, FALSE, msec, QS_ALLEVENTS );
#endif
}


/* -- INTERFACE FUNCTIONS -------------------------------------------------- */


/*
 * Enters the freeglut processing loop.
 * Stays until the "ExecState" changes to "GLUT_EXEC_STATE_STOP".
 */
void FGAPIENTRY glutMainLoopMT( void )
{
    //int action; //MT

#if TARGET_HOST_WIN32 || TARGET_HOST_WINCE
    SFG_Window *window = (SFG_Window *)fgStructure.Windows.First ;
#endif

    FREEGLUT_EXIT_IF_NOT_INITIALISED ( "glutMainLoop" );

#if TARGET_HOST_WIN32 || TARGET_HOST_WINCE
    /*
     * Processing before the main loop:  If there is a window which is open and
     * which has a visibility callback, call it.  I know this is an ugly hack,
     * but I'm not sure what else to do about it.  Ideally we should leave
     * something uninitialized in the create window code and initialize it in
     * the main loop, and have that initialization create a "WM_ACTIVATE"
     * message.  Then we would put the visibility callback code in the
     * "case WM_ACTIVATE" block below.         - John Fay -- 10/24/02
     */
    while( window )
    {
        if ( FETCH_WCB( *window, Visibility ) )
        {
            SFG_Window *current_window = fgStructure.CurrentWindow ;

            INVOKE_WCB( *window, Visibility, ( window->State.Visible ) );
            fgSetWindow( current_window );
        }

        window = (SFG_Window *)window->Node.Next ;
    }
#endif

    fgState.ExecState = GLUT_EXEC_STATE_RUNNING ;
    while( fgState.ExecState == GLUT_EXEC_STATE_RUNNING )
    {
        SFG_Window *window;

        glutMainLoopEvent( );
        /*
         * Step through the list of windows, seeing if there are any
         * that are not menus
         */
        for( window = ( SFG_Window * )fgStructure.Windows.First;
             window;
             window = ( SFG_Window * )window->Node.Next )
            if ( ! ( window->IsMenu ) )
                break;

        if( ! window )
            fgState.ExecState = GLUT_EXEC_STATE_STOP;
        else
        {
            if( fgState.IdleCallback )
            {
                if( fgStructure.CurrentWindow &&
                    fgStructure.CurrentWindow->IsMenu )
                    /* fail safe */
                    fgSetWindow( window );
                fgState.IdleCallback( );
            }

            fghSleepForEvents( );
        }
    }

#if 0 //Marc Toussaint
    /*
     * When this loop terminates, destroy the display, state and structure
     * of a freeglut session, so that another glutInit() call can happen
     *
     * Save the "ActionOnWindowClose" because "fgDeinitialize" resets it.
     */
    action = fgState.ActionOnWindowClose;
    fgDeinitialize( );
    if( action == GLUT_ACTION_EXIT )
        exit( 0 );
#endif
}




/*** END OF FILE ***/
