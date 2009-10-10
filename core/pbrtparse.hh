/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     STRING = 258,
     ID = 259,
     NUM = 260,
     LBRACK = 261,
     RBRACK = 262,
     ACCELERATOR = 263,
     ACTIVETRANSFORM = 264,
     ALL = 265,
     AREALIGHTSOURCE = 266,
     ATTRIBUTEBEGIN = 267,
     ATTRIBUTEEND = 268,
     CAMERA = 269,
     CONCATTRANSFORM = 270,
     COORDINATESYSTEM = 271,
     COORDSYSTRANSFORM = 272,
     ENDTIME = 273,
     FILM = 274,
     IDENTITY = 275,
     LIGHTSOURCE = 276,
     LOOKAT = 277,
     MAKENAMEDMATERIAL = 278,
     MATERIAL = 279,
     NAMEDMATERIAL = 280,
     OBJECTBEGIN = 281,
     OBJECTEND = 282,
     OBJECTINSTANCE = 283,
     PIXELFILTER = 284,
     RENDERER = 285,
     REVERSEORIENTATION = 286,
     ROTATE = 287,
     SAMPLER = 288,
     SCALE = 289,
     SHAPE = 290,
     STARTTIME = 291,
     SURFACEINTEGRATOR = 292,
     TEXTURE = 293,
     TRANSFORMBEGIN = 294,
     TRANSFORMEND = 295,
     TRANSFORMTIMES = 296,
     TRANSFORM = 297,
     TRANSLATE = 298,
     VOLUME = 299,
     VOLUMEINTEGRATOR = 300,
     WORLDBEGIN = 301,
     WORLDEND = 302,
     HIGH_PRECEDENCE = 303
   };
#endif
/* Tokens.  */
#define STRING 258
#define ID 259
#define NUM 260
#define LBRACK 261
#define RBRACK 262
#define ACCELERATOR 263
#define ACTIVETRANSFORM 264
#define ALL 265
#define AREALIGHTSOURCE 266
#define ATTRIBUTEBEGIN 267
#define ATTRIBUTEEND 268
#define CAMERA 269
#define CONCATTRANSFORM 270
#define COORDINATESYSTEM 271
#define COORDSYSTRANSFORM 272
#define ENDTIME 273
#define FILM 274
#define IDENTITY 275
#define LIGHTSOURCE 276
#define LOOKAT 277
#define MAKENAMEDMATERIAL 278
#define MATERIAL 279
#define NAMEDMATERIAL 280
#define OBJECTBEGIN 281
#define OBJECTEND 282
#define OBJECTINSTANCE 283
#define PIXELFILTER 284
#define RENDERER 285
#define REVERSEORIENTATION 286
#define ROTATE 287
#define SAMPLER 288
#define SCALE 289
#define SHAPE 290
#define STARTTIME 291
#define SURFACEINTEGRATOR 292
#define TEXTURE 293
#define TRANSFORMBEGIN 294
#define TRANSFORMEND 295
#define TRANSFORMTIMES 296
#define TRANSFORM 297
#define TRANSLATE 298
#define VOLUME 299
#define VOLUMEINTEGRATOR 300
#define WORLDBEGIN 301
#define WORLDEND 302
#define HIGH_PRECEDENCE 303




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 147 "core/pbrtparse.yy"
{
char string[1024];
float num;
ParamArray *ribarray;
}
/* Line 1529 of yacc.c.  */
#line 151 "core/pbrtparse.hpp"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

