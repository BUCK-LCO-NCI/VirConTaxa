/* Copyright (C) 2000 MySQL AB

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

/*
  Code for handling strings with can grow dynamicly.
  Copyright Monty Program KB.
  By monty.
*/
#include <sys/types.h>
#include <string.h>
#include <stdlib.h>


#define CALLER_INFO_PROTO   , const char *sFile, uint uLine
#define CALLER_INFO         , __FILE__, __LINE__
#define ORIG_CALLER_INFO    , sFile, uLine
#define TRUE 1
#define FALSE 0
#define MY_WME  16
#define DBUG_RETURN(a1) return(a1)
#define DBUG_ENTER(a1)

typedef char my_bool;
typedef unsigned long ulong;

typedef struct st_dynamic_string
{
  char *str;
  uint length,max_length,alloc_increment;
} DYNAMIC_STRING;



my_bool dynstr_append(DYNAMIC_STRING *str, const char *append);
my_bool dynstr_append_mem(DYNAMIC_STRING *str, const char *append,
			  uint length);
my_bool dynstr_set(DYNAMIC_STRING *str, const char *init_str);
my_bool dynstr_realloc(DYNAMIC_STRING *str, ulong additional_size);
void dynstr_free(DYNAMIC_STRING *str);


my_bool init_dynamic_string(DYNAMIC_STRING *str, const char *init_str,
			    uint init_alloc, uint alloc_increment)
{
  uint length;
  DBUG_ENTER("init_dynamic_string");

  if (!alloc_increment)
    alloc_increment=128;
  length=1;
  if (init_str && (length= (uint) strlen(init_str)+1) < init_alloc)
    init_alloc=((length+alloc_increment-1)/alloc_increment)*alloc_increment;
  if (!init_alloc)
    init_alloc=alloc_increment;

  if (!(str->str=(char*) malloc(init_alloc)))
    DBUG_RETURN(TRUE);
  str->length=length-1;
  if (init_str)
    memcpy(str->str,init_str,length);
  str->max_length=init_alloc;
  str->alloc_increment=alloc_increment;
  DBUG_RETURN(FALSE);
}


my_bool dynstr_set(DYNAMIC_STRING *str, const char *init_str)
{
  uint length=0;
  DBUG_ENTER("dynstr_set");

  if (init_str && (length= (uint) strlen(init_str)+1) > str->max_length)
  {
    str->max_length=((length+str->alloc_increment-1)/str->alloc_increment)*
      str->alloc_increment;
    if (!str->max_length)
      str->max_length=str->alloc_increment;
    if (!(str->str=(char*) realloc(str->str,str->max_length)))
      DBUG_RETURN(TRUE);
  }
  if (init_str)
  {
    str->length=length-1;
    memcpy(str->str,init_str,length);
  }
  else
    str->length=0;
  DBUG_RETURN(FALSE);
}


my_bool dynstr_realloc(DYNAMIC_STRING *str, ulong additional_size)
{
  DBUG_ENTER("dynstr_realloc");

  if (!additional_size) DBUG_RETURN(FALSE);
  if (str->length + additional_size > str->max_length)
  {
    str->max_length=((str->length + additional_size+str->alloc_increment-1)/
		     str->alloc_increment)*str->alloc_increment;
    if (!(str->str=(char*) realloc(str->str,str->max_length)))
      DBUG_RETURN(TRUE);
  }
  DBUG_RETURN(FALSE);
}


my_bool dynstr_append(DYNAMIC_STRING *str, const char *append)
{
  return dynstr_append_mem(str,append,strlen(append));
}


my_bool dynstr_append_mem(DYNAMIC_STRING *str, const char *append,
			  uint length)
{
  char *new_ptr;
  if (str->length+length >= str->max_length)
  {
    uint new_length=(str->length+length+str->alloc_increment)/
      str->alloc_increment;
    new_length*=str->alloc_increment;
    if (!(new_ptr=(char*) realloc(str->str,new_length)))
      return TRUE;
    str->str=new_ptr;
    str->max_length=new_length;
  }
  memcpy(str->str + str->length,append,length);
  str->length+=length;
  str->str[str->length]=0;			/* Safety for C programs */
  return FALSE;
}


void dynstr_free(DYNAMIC_STRING *str)
{
  if (str->str)
  {
    free(str->str);
    str->str=0;
  }
}
