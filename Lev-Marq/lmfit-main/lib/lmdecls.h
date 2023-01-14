/*
 * Library:   lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:      decls.h
 *
 * Contents:  Macros for C in C++ and for Windows import/export
 *
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2004-2013)
 *
 * License:   see ../COPYING (FreeBSD)
 *
 * Homepage:  https://jugit.fz-juelich.de/mlz/lmfit
 */

#ifndef LM_DLL /* will be defined below */

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS /* empty */
#define __END_DECLS   /* empty */
#endif /* __cplusplus */

#if WIN32
#ifdef LMFIT_EXPORT
#define LM_DLL __declspec(dllexport)
#else
#define LM_DLL __declspec(dllimport)
#endif
#else
#define LM_DLL
#endif /* WIN32 */

#endif /* LM_DLL */
