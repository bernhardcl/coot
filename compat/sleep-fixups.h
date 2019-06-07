/* compat/sleep-fixups.h
 * 
 * Copyright 2009 by Bernhard Lohkamp
 * Author: Bernhard Lohkamp
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#if defined (WINDOWS_MINGW) || defined (_MSC_VER)
#ifdef sleep
#sleep
#endif
#ifdef usleep
#usleep
#endif
#if !defined (usleep) || !defined (sleep)
#ifndef Sleep
#include <windows.h>
#endif /* Sleep */
#ifndef usleep
#define usleep(t) Sleep(t/1000)
#endif /* usleep */
#ifndef sleep
#define sleep(t) Sleep(1000*t)
#endif /* sleep */
#endif /* usleep */
#endif /* windows */

