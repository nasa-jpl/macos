/* giza - a scientific plotting library built on cairo
 *
 * Copyright (c) 2010      James Wetter and Daniel Price
 * Copyright (c) 2010-2015 Daniel Price
 *
 * This library is free software; and you are welcome to redistribute
 * it under the terms of the GNU General Public License
 * (GPL, see LICENSE file for details) and the provision that
 * this notice remains intact. If you modify this file, please
 * note section 2a) of the GPLv2 states that:
 *
 *  a) You must cause the modified files to carry prominent notices
 *     stating that you changed the files and the date of any change.
 *
 * This software is distributed "AS IS", with ABSOLUTELY NO WARRANTY.
 * See the GPL for specific language governing rights and limitations.
 *
 * The Original code is the giza plotting library.
 *
 * Contributor(s):
 *      James Wetter <wetter.j@gmail.com>
 *      Daniel Price <daniel.price@monash.edu> (main contact)
 */

#include "giza-private.h"
#include "giza-drivers-private.h"
#include "giza-driver-xw-private.h"
#include "giza-driver-win-private.h"
#include "giza-driver-png-private.h"
#include "giza-driver-pdf-private.h"
#include "giza-driver-ps-private.h"
#include "giza-driver-eps-private.h"
#include "giza-driver-svg-private.h"
#include "giza-driver-null-private.h"
#include "giza-io-private.h"
#include "giza-viewport-private.h"
#include "giza-window-private.h"
#include "giza-transforms-private.h"
#include "giza-character-size-private.h"
#include "giza-text-private.h"
#include "giza-set-font-private.h"
#include "giza-colour-private.h"
#include "giza-arrow-style-private.h"
#include "giza-fill-private.h"
#include "giza-line-style-private.h"
#include "giza-band-private.h"
#include "giza-stroke-private.h"
#include "giza-subpanel-private.h"
#include "giza-text-background-private.h"
#include <giza.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>

#define GIZA_DEFAULT_MARGIN 0

static void _giza_set_prefix (const char *prefix);
static int _giza_get_internal_id(int devid);
static void _giza_close_device_unchecked (void);
static void _giza_init_device_struct(giza_device_t*);

/* global settings */
giza_settings_t Sets;
int id;
/* per-device settings */
giza_device_t Dev[GIZA_MAX_DEVICES];

int
giza_open_device (const char *newDeviceName, const char *newPrefix)
{
  return giza_open_device_size (newDeviceName, newPrefix, 0., 0., 0);
}

int
giza_open_device_size (const char *newDeviceName, const char *newPrefix, double width, double height, int units)
{
  static int didInit = 0;

  if( !didInit ) {
	  giza_device_t* pDev;
      for(pDev = &Dev[0]; pDev < &Dev[GIZA_MAX_DEVICES]; pDev++)
          _giza_init_device_struct( pDev );
       didInit = 1;
  }

  if( !newDeviceName || !strlen(newDeviceName) ) {
      _giza_error("giza_open_device", newDeviceName==NULL ? "(nullptr devicename)" : "empty device name");
      return -1;
  }

  int   newId;
  for( newId=0; newId<GIZA_MAX_DEVICES; newId++ )
      if( Dev[newId].deviceOpen==0 )
          break;
  if( newId>=GIZA_MAX_DEVICES ) {
      _giza_error("giza_open_device", "No more free devices (%d in use). Close some first.", GIZA_MAX_DEVICES);
      return -1;
  }

  const int prevId = id;
  id = newId;
  _giza_init();
  _giza_set_prefix( newPrefix!=0 ? newPrefix : GIZA_DEFAULT_PREFIX );
  Dev[id].number_format = GIZA_NUMBER_FORMAT_AUTO;

  char firstchar = newDeviceName[0];
  if (firstchar == '?')
    {
    Dev[id].type = _giza_prompt_for_device ();
    }
  else
    {
      char const *devTypeStr;
      _giza_split_device_string (newDeviceName, &devTypeStr);
      Dev[id].type = _giza_device_to_int (devTypeStr);
    }

  int success;
  switch (Dev[id].type)
    {
#ifdef _GIZA_HAS_XW
    case GIZA_DEVICE_XW:
      success = _giza_open_device_xw (width, height, units);
      break;
#endif
#ifdef _GIZA_HAS_WIN
    case GIZA_DEVICE_WIN:
      success = _giza_open_device_win (width, height, units);
      break;
#endif
    case GIZA_DEVICE_PNG:
      success = _giza_open_device_png (width, height, units);
      break;
    case GIZA_DEVICE_PDF:
      success = _giza_open_device_pdf (width, height, units, 0);
      break;
    case GIZA_DEVICE_VPDF:
      success = _giza_open_device_pdf(width, height, units, 1);
      break;
    case GIZA_DEVICE_PS:
      success = _giza_open_device_ps (width, height, units, 0);
      break;
    case GIZA_DEVICE_VPS:
      success = _giza_open_device_ps (width, height, units, 1);
      break;
    case GIZA_DEVICE_SVG:
      success = _giza_open_device_svg (width, height, units, 0);
      break;
    case GIZA_DEVICE_NULL:
      success = _giza_open_device_null (width, height, units);
      break;
#ifdef _GIZA_HAS_EPS
    case GIZA_DEVICE_EPS:
      success = _giza_open_device_eps (width, height, units, 0);
      break;
#endif
    case GIZA_DEVICE_IV:
    default:
      _giza_error ("giza_open_device", "Unknown device, device not opened");
      success = 666;
      break;
    }

  if (success!=0 || Dev[id].surface==0) {
      _giza_close_device_unchecked ();
      id = prevId;
      return -1;
  }
  Dev[id].context = cairo_create (Dev[id].surface);
  if (!Dev[id].context)
    {
      _giza_error ("giza_open_device", "Could not create cairo context.");
      _giza_close_device_unchecked();
      id = prevId;
      return -1;
    }

  Dev[id].deviceOpen = 1;
  Dev[id].defaultBackgroundAlpha = 1.;

  giza_set_text_background (-1);
  giza_start_prompting ();
  _giza_init_subpanel ();
  _giza_init_arrow_style ();
  _giza_init_line_style ();
  _giza_init_colour_index ();

  giza_draw_background ();
  _giza_init_colour_table ();
  _giza_set_trans (GIZA_TRANS_IDEN);
  _giza_init_norm ();

  _giza_init_character_height ();
  _giza_init_font ();
  _giza_init_text_background();

  _giza_init_window ();
  giza_set_viewport_default ();
  giza_set_line_width (1);

  _giza_init_fill ();
  _giza_init_band_style ();
  _giza_init_save ();
  giza_set_clipping(1);
  return id+1;
}

int
giza_open_device_size_float (const char *newDeviceName, const char *newPrefix, float width, float height, int units)
{
  return giza_open_device_size (newDeviceName, newPrefix, (double) width, (double) height, units);
}

int _giza_get_internal_id(int devid)
{
  return devid - 1;
}

void
giza_get_device_id (int *devid)
{
  *devid = id + 1;
  return;
}

void
giza_select_device (int devid)
{
    int tmpid = _giza_get_internal_id(devid);
    if( tmpid<0 || tmpid>=GIZA_MAX_DEVICES ) {
        _giza_error ("giza_select_device", "Invalid device number %d selected", devid);
        return;
    }
    if (!Dev[tmpid].deviceOpen) {
        _giza_error ("giza_select_device", "Invalid/closed device %d selected", devid);
        return;
    }
    id = tmpid;
}

void
giza_flush_device (void)
{
  if (!_giza_check_device_ready ("giza_flush_device"))
    return;

  Dev[id].drawn = 1;

  if (!Dev[id].buf)
    {
     switch (Dev[id].type)
       {
#ifdef _GIZA_HAS_XW
       case GIZA_DEVICE_XW:
         _giza_flush_device_xw ();
         break;
#endif
#ifdef _GIZA_HAS_WIN
       case GIZA_DEVICE_WIN:
         _giza_flush_device_win ();
         break;
#endif
       default:
         if (!Dev[id].surface)
           {
             _giza_error ("giza_flush_device", "No device open, cannot flush");
             return;
           } else {
             cairo_surface_flush(Dev[id].surface);
           }
         return;
       }
    }
  return;
}

void
giza_change_page (void)
{
  if (Dev[id].nx > 1 || Dev[id].ny > 1) {
     int newpage;
     _giza_advance_panel(&newpage);
     if (!newpage) return;
  }

  if (!Dev[id].drawn && !Dev[id].resize) {
    giza_draw_background ();
    return;
  }

  if (!_giza_check_device_ready ("giza_change_page"))
    return;

  if (Dev[id].prompting && Dev[id].isInteractive && !Dev[id].firstPage)
    {
      _giza_newpage_prompt();
    }
  Dev[id].firstPage = 0;
  giza_save();

  switch (Dev[id].type)
    {
#ifdef _GIZA_HAS_XW
    case GIZA_DEVICE_XW:
      _giza_change_page_xw ();
      break;
#endif
#ifdef _GIZA_HAS_WIN
    case GIZA_DEVICE_WIN:
      _giza_change_page_win ();
      break;
#endif
    case GIZA_DEVICE_PNG:
      _giza_change_page_png ();
      break;
    case GIZA_DEVICE_SVG:
      _giza_change_page_svg ();
      break;
    case GIZA_DEVICE_PDF:
    case GIZA_DEVICE_VPDF:
      _giza_change_page_pdf ();
      break;
    case GIZA_DEVICE_PS:
    case GIZA_DEVICE_VPS:
      _giza_change_page_ps ();
      break;
    case GIZA_DEVICE_NULL:
      _giza_change_page_null ();
      break;
#ifdef _GIZA_HAS_EPS
    case GIZA_DEVICE_EPS:
      _giza_change_page_eps ();
      break;
#endif
    default:
      _giza_error ("giza_change_page", "No device open");
      return;
    }

  if (Dev[id].drawn)
      Dev[id].pgNum++;

  giza_set_panel(1,1);
  giza_set_window (Dev[id].Win.xmin, Dev[id].Win.xmax,
                   Dev[id].Win.ymin, Dev[id].Win.ymax);
  giza_restore();

  Dev[id].drawn = 0;
  Dev[id].resize = 0;
  giza_draw_background ();
  giza_flush_device ();
  return;
}

void
giza_close_devices (void)
{
  for(id = 0; id<GIZA_MAX_DEVICES; id++)
      if( Dev[id].deviceOpen )
          giza_close_device ();
  _giza_unload_font_cache ();
}

void
giza_close_device (void)
{
  if (!_giza_check_device_ready ("giza_close_device"))
    return;

  if (Dev[id].prompting && Dev[id].isInteractive)
    {
      cairo_surface_finish (Dev[id].surface);
      _giza_newpage_prompt();
    }
  _giza_close_device_unchecked();

  unsigned int nOpen = 0;
  giza_device_t* p;
  for(p=&Dev[0]; p<&Dev[GIZA_MAX_DEVICES]; p++)
      if( p->deviceOpen )
          nOpen++;
  if( nOpen==0 )
    _giza_unload_font_cache ();
}

void
_giza_close_device_unchecked (void) {
  if (Dev[id].context)
      cairo_destroy(Dev[id].context);

  switch (Dev[id].type)
    {
#ifdef _GIZA_HAS_XW
    case GIZA_DEVICE_XW:
      _giza_close_device_xw ();
      break;
#endif
#ifdef _GIZA_HAS_WIN
    case GIZA_DEVICE_WIN:
      _giza_close_device_win ();
      break;
#endif
    case GIZA_DEVICE_PNG:
      _giza_close_device_png (1);
      break;
    case GIZA_DEVICE_SVG:
      _giza_close_device_svg ();
      break;
    case GIZA_DEVICE_PDF:
    case GIZA_DEVICE_VPDF:
      _giza_close_device_pdf ();
      break;
    case GIZA_DEVICE_PS:
    case GIZA_DEVICE_VPS:
      _giza_close_device_ps ();
      break;
    case GIZA_DEVICE_NULL:
      _giza_close_device_null ();
      break;
#ifdef _GIZA_HAS_EPS
    case GIZA_DEVICE_EPS:
      _giza_close_device_eps ();
      break;
#endif
    default:
       break;
    }

  _giza_free_font ();
  _giza_free_colour_table ();
  _giza_init_device_struct( &Dev[id] );
}

void
_giza_init_device_struct(giza_device_t* ptrDev) {
    memset(ptrDev, 0x0, sizeof(giza_device_t));
    ptrDev->type = GIZA_DEVICE_IV;
}

int
giza_query_device (const char *querytype, char *returnval, int* rlen)
{
  int       ierr = 0;
  char      devType[12];
  const int max_chars = *rlen - 1;

  if (max_chars<=0)
    {
        _giza_warning("giza_query_device", "destination string says it has %d characters available querying %s", max_chars, querytype);
        return 1;
    }

  if (!strcasecmp(querytype,"type"))
     {
       if (!_giza_int_to_device(Dev[id].type,returnval, max_chars))
         {
           ierr = 1;
         }
     }
  else if (!strcasecmp(querytype,"cursor"))
    {
      strncpy(returnval, Dev[id].isInteractive ? "YES" : "NO", max_chars);
    }
  else if (!strcasecmp(querytype,"hardcopy"))
    {
      strncpy(returnval, Dev[id].isInteractive ? "NO" : "YES", max_chars);
    }
  else if (!strcasecmp(querytype,"state"))
    {
      strncpy(returnval, Dev[id].deviceOpen ? "OPEN" : "CLOSED", max_chars);
    }
  else if (!strcasecmp(querytype,"user"))
    {
       strncpy(returnval,getlogin(),max_chars);
    }
  else if (!strcasecmp(querytype,"device"))
    {
       strncpy(returnval,Dev[id].prefix,max_chars);
    }
  else if (!strcasecmp(querytype,"dev/type"))
    {
       strncpy(returnval,Dev[id].prefix,max_chars);
       if (!_giza_int_to_device(Dev[id].type,devType,sizeof(devType)))
         {
           ierr = 1;
           strncat(returnval,devType,4*sizeof(char));
         }
    }
  else if (!strcasecmp(querytype,"file"))
    {
       if (!Dev[id].isInteractive)
         {
           strncpy(returnval,Dev[id].prefix,max_chars);
         }
       else
         {
           strncpy(returnval, " ", max_chars);
         }
    }
  returnval[max_chars] = '\0';
  return ierr;
}

int
_giza_get_key_press (int mode, int moveCurs, int nanc, const double *xanch, const double *yanch, double *x, double *y, char *ch)
{
  if (!_giza_check_device_ready ("_giza_get_key_press"))
    return 1;

  switch (Dev[id].type)
    {
    case GIZA_DEVICE_NULL:
    case GIZA_DEVICE_PDF:
    case GIZA_DEVICE_VPDF:
    case GIZA_DEVICE_PNG:
    case GIZA_DEVICE_SVG:
    case GIZA_DEVICE_PS:
    case GIZA_DEVICE_VPS:
      _giza_warning ("giza_get_key_press", "Current device does not support a cursor, returning x = 0, y = 0, ch = a");
      return 1;
      break;
#ifdef _GIZA_HAS_XW
    case GIZA_DEVICE_XW:
      _giza_get_key_press_xw (mode, moveCurs, nanc, xanch, yanch, x, y, ch);
      return 0;
      break;
#endif
#ifdef _GIZA_HAS_WIN
    case GIZA_DEVICE_WIN:
      _giza_get_key_press_win (mode, moveCurs, nanc, xanch, yanch, x, y, ch);
      return 0;
      break;
#endif
    case GIZA_DEVICE_IV:
    default:
      _giza_error ("giza_get_key_press", "Unknown device");
      return 1;
    }
}

void
_giza_split_device_string (const char *deviceString, char const **devType)
{
  if (deviceString == NULL)
    return;

  *devType = strrchr (deviceString, (int) '/');

  if (*devType == NULL)
    {
      *devType = strrchr (deviceString, (int) '.');
      if (*devType == NULL)
         {
           *devType = deviceString;
           return;
         }
    }

  int lendev  = strlen(deviceString);
  int lentype = strlen(*devType);
  int nameLength = lendev - lentype;
  if (nameLength != 0)
    {
      char devName[nameLength + 1];
      strncpy (devName, deviceString, (size_t) nameLength);
      devName[nameLength] = '\0';
      _giza_set_prefix (devName);
    }
}

int
_giza_device_to_int (const char *newDeviceName)
{
  int newDevice;

  char devName[strlen (newDeviceName) + 1];

  int i;
  for (i = 0; i < (signed ) strlen (newDeviceName); ++i)
    {
      devName[i] = (char) tolower (newDeviceName[i]);
    }
  devName[i] = '\0';

  if (!strcmp (devName, "/null"))
    newDevice = GIZA_DEVICE_NULL;
#ifdef _GIZA_HAS_XW
  else if (!strcmp (devName, "/xw")
        || !strcmp (devName, "/xs")
        || !strcmp (devName, "/xserve")
        || !strcmp (devName, "/xwindow")
        || !strcmp (devName, "/xwin"))
    newDevice = GIZA_DEVICE_XW;
#endif
#ifdef _GIZA_HAS_WIN
  else if (!strcmp (devName, "/xw")
        || !strcmp (devName, "/xs")
        || !strcmp (devName, "/xserve")
        || !strcmp (devName, "/xwindow")
        || !strcmp (devName, "/xwin")
        || !strcmp (devName, "/win"))
    newDevice = GIZA_DEVICE_WIN;
#endif
  else if (!strcmp (devName, "/png"))
    newDevice = GIZA_DEVICE_PNG;
  else if (!strcmp (devName, "/svg"))
    newDevice = GIZA_DEVICE_SVG;
  else if (!strcmp (devName, "/pdf"))
    newDevice = GIZA_DEVICE_PDF;
  else if (!strcmp (devName, "/ps")
        || !strcmp (devName, "/cps")
        || !strcmp (devName, "/postscript"))
    newDevice = GIZA_DEVICE_PS;
  else if (!strcmp (devName, "/vpdf"))
    newDevice = GIZA_DEVICE_VPDF;
  else if (!strcmp (devName, "/vps"))
    newDevice = GIZA_DEVICE_VPS;
#ifdef _GIZA_HAS_EPS
  else if (!strcmp (devName, "/eps"))
    newDevice = GIZA_DEVICE_EPS;
#endif
  else
    {
      if (!strcmp (devName, ".png"))
        newDevice = GIZA_DEVICE_PNG;
      else if (!strcmp (devName, ".svg"))
        newDevice = GIZA_DEVICE_SVG;
      else if (!strcmp (devName, ".pdf"))
        newDevice = GIZA_DEVICE_PDF;
      else if (!strcmp (devName, ".ps"))
        newDevice = GIZA_DEVICE_PS;
#ifdef _GIZA_HAS_EPS
      else if (!strcmp (devName, ".eps"))
        newDevice = GIZA_DEVICE_EPS;
#endif
      else {
         newDevice = GIZA_DEVICE_IV;
      }
    }

  return newDevice;
}

int
_giza_int_to_device (int numDevice, char *DeviceName, int rval)
{
 int       ierr = 0;
 const int max_chars = rval - 1;
 switch (numDevice)
    {
    case GIZA_DEVICE_NULL:
      strncpy(DeviceName,"/null",max_chars);
      break;
    case GIZA_DEVICE_PDF:
      strncpy(DeviceName,"/pdf",max_chars);
      break;
    case GIZA_DEVICE_VPDF:
      strncpy(DeviceName,"/vpdf",max_chars);
      break;
    case GIZA_DEVICE_PNG:
      strncpy(DeviceName,"/png",max_chars);
      break;
    case GIZA_DEVICE_SVG:
      strncpy(DeviceName,"/svg",max_chars);
      break;
    case GIZA_DEVICE_PS:
      strncpy(DeviceName,"/ps",max_chars);
      break;
    case GIZA_DEVICE_VPS:
      strncpy(DeviceName,"/vps",max_chars);
      break;
#ifdef _GIZA_HAS_XW
    case GIZA_DEVICE_XW:
      strncpy(DeviceName,"/xw",max_chars);
      break;
#endif
#ifdef _GIZA_HAS_WIN
    case GIZA_DEVICE_WIN:
      strncpy(DeviceName,"/xw",max_chars);
      break;
#endif
#ifdef _GIZA_HAS_EPS
    case GIZA_DEVICE_EPS:
      strncpy(DeviceName,"/eps",max_chars);
      break;
#endif
    default:
      strncpy(DeviceName," ",max_chars);
      ierr = 1;
    }
    DeviceName[max_chars] = '\0';
    return ierr;
}

void
giza_print_device_list (void)
{
   _giza_display_devices();
}

void
_giza_init_device_list (char **deviceList)
{
  *deviceList = malloc (1000 * sizeof(char));
  *deviceList[0] = '\0';
#ifdef _GIZA_HAS_XW
  strcat (*deviceList, "Interactive devices:\n\t/xw\t(X Window)\n");
#endif
#ifdef _GIZA_HAS_WIN
  strcat (*deviceList, "Interactive devices:\n\t/xw\t(Win32 Window)\n");
#endif
  strcat (*deviceList, "Non-interactive file formats:\n");
  strcat (*deviceList, "\t/png or file.png   (Portable network graphics file)\n");
  strcat (*deviceList, "\t/svg or file.svg   (Scalable vector graphics file)\n");
#ifdef _GIZA_HAS_EPS
  strcat (*deviceList, "\t/eps or file.eps   (Encapsulated Postscript, one file per page)\n");
#endif
  strcat (*deviceList, "\t/pdf or file.pdf   (Portable document format file)\n");
  strcat (*deviceList, "\t/vpdf              (Portable document format file portrait)\n");
  strcat (*deviceList, "\t/ps  or file.ps    (Postscript file, multiple pages per file)\n");
  strcat (*deviceList, "\t/vps               (Postscript file portrait, multiple pages per file)\n");
  strcat (*deviceList, "\t/null              (Null device)\n");
}

void
_giza_free_device_list (char *deviceList)
{
  if (deviceList)
    {
      free (deviceList);
    }
}

void
_giza_init_norm (void)
{
  if (!_giza_check_device_ready ("_giza_init_norm"))
    return;

  double xx,yy,x0,y0;

  switch (Dev[id].type)
    {
#ifdef _GIZA_HAS_XW
    case GIZA_DEVICE_XW:
      _giza_init_norm_xw ();
      break;
#endif
#ifdef _GIZA_HAS_WIN
    case GIZA_DEVICE_WIN:
      /* Win device uses same margin-based norm as default — fall through */
#endif
    default:
      xx = (double) Dev[id].width - 2.*GIZA_DEFAULT_MARGIN;
      yy = 2.*GIZA_DEFAULT_MARGIN - (double) Dev[id].height;
      x0 = GIZA_DEFAULT_MARGIN;
      y0 = (double) Dev[id].height - GIZA_DEFAULT_MARGIN;
      cairo_matrix_init (&(Dev[id].Win.normCoords), xx, 0, 0, yy, x0, y0);
      _giza_set_trans (GIZA_TRANS_NORM);
      break;
    }
}

void
_giza_expand_clipping (void)
{
    if (!_giza_check_device_ready ("_giza_expand_clipping"))
    return;
  switch (Dev[id].type)
    {
#ifdef _GIZA_HAS_XW
    case GIZA_DEVICE_XW:
      _giza_expand_clipping_xw ();
      break;
#endif
#ifdef _GIZA_HAS_WIN
    case GIZA_DEVICE_WIN:
      _giza_expand_clipping_win ();
      break;
#endif
    default:
      _giza_set_trans (GIZA_TRANS_IDEN);
      cairo_reset_clip (Dev[id].context);
      cairo_rectangle (Dev[id].context, 0, 0, Dev[id].width, Dev[id].height);
      cairo_clip (Dev[id].context);
      break;
    }
}

void
_giza_set_prefix (const char *prefix)
{
  if (strlen(prefix) > sizeof(Dev[id].prefix))
     {
        _giza_error("giza_set_prefix","device name exceeds maximum string length");
     }
  strncpy (Dev[id].prefix, prefix, sizeof(Dev[id].prefix));
}

int
_giza_init_band (int mode)
{
  if (mode == 0) return 0;

  int success = 1;
  switch (Dev[id].type)
    {
#ifdef _GIZA_HAS_XW
      case GIZA_DEVICE_XW:
        success = _giza_init_band_xw ();
        break;
#endif
#ifdef _GIZA_HAS_WIN
      case GIZA_DEVICE_WIN:
        success = _giza_init_band_win ();
        break;
#endif
      default:
        _giza_error ("_giza_init_band", "band not implemented for this device");
	break;
    }
  _giza_set_line_style (Band.ls, Band.box);
  double lwDevice = Band.lw * Dev[id].deviceUnitsPermm * 0.25;
  cairo_set_line_width (Band.box, lwDevice);

  return success;
}

void _giza_lowercase(const char *string, char *lowerstring)
{
   int  i = 0;
   while ( string[i] )
   {
      lowerstring[i] = (char) tolower(string[i]);
      i++;
   }
   lowerstring[i] = '\0';
   return;
}

void _giza_trim(char *str) {
   int len = strlen(str);
   int i;
   for (i=0; i<len && isspace(str[i]); i++, str++);
   for (i=len-1; i>=0 && isspace(str[i]); str[i]=0, i--);
   return;
}

void _giza_get_filename_for_device (char *filename, char *prefix, int pgNum, char *extension,
                                    int lastpage)
{
  char lprefix[strlen(prefix) + 1];
  char lextens[strlen(extension) + 1];
  _giza_lowercase(prefix,lprefix);
  _giza_lowercase(extension,lextens);

  char *prefixtrim = prefix;
  char *ext = extension;
  _giza_trim(prefixtrim);
  _giza_trim(ext);

  if (!strstr(lprefix,lextens)) {
     if (pgNum == 0 && lastpage != 0) {
        sprintf (filename, "%s%s", prefixtrim, ext);
     } else {
        sprintf (filename, "%s_%04d%s", prefixtrim, pgNum, ext);
     }
  } else {
     if (pgNum == 0 && lastpage != 0) {
        sprintf (filename, "%s", prefixtrim);
     } else {
        char *firstpart = strsep (&prefixtrim,".");
        sprintf (filename, "%s_%04d%s", firstpart, pgNum, ext);
     }
  }
}
